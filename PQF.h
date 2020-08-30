//
//  PQF.h
//  RS Quotient Filter
//

#ifndef PQF_h
#define PQF_h




class PQF {

    static constexpr int batch_size = 128;
    
public:
    unsigned char * ar;  // the raw data
    int    num_lines; // The number of 64 byte blocks.
    int log_num_lines; // the number of bits in indexing blocks
    unsigned long rm;
    unsigned long long int  hs[256];
    unsigned long long int  hs2[256];
    unsigned long long int  hs3[256];
    
    PQF(int sz){
        if (posix_memalign((void **)&ar, 64, sz * 64  ) != 0) {
            cout << "error on memalign" << endl;
            exit(1);
        }
        bzero(ar,  sz * 64);
        num_lines = sz;
        log_num_lines = log2(sz);
        rm = 0xffffffffffff; // A mask of 48 bits.
        // hs is an involution on log_num_lines + 8 bits.
        // We
        for (int i =0;i < 256; i++){
            hs[i] = ( (unsigned  int)hash1(i)&( (1UL<<(log_num_lines+8))-1  ));
            while ((hs[i]&(num_lines-1)) == 0)    hs[i] = ( (unsigned  int)hash1( rand() +1 )&( (1UL<<(log_num_lines+8))-1  ))  ;
        }
        // hs2 is for lines
        for (int i =0;i < 256; i++){
            hs2[i] = ( (unsigned  int)hash1( rand() +1 )&(num_lines- 1  ));
            while ((hs2[i]&(num_lines - 1 )) == 0 )
                hs2[i] = ( (unsigned  int)hash1( rand() +1 )&(  num_lines - 1))   ;
        }
        
        //hs3 is for places (0--127)
        for (int i =0;i < 256; i++){
            hs3[i] = ( (unsigned  int)hash1(rand() )&( 63  ));
            while ((hs3[i]&(63))  == 0  )
                hs3[i] = ( (unsigned  int)hash1( rand() +1 )&( 63))   ;
        }
        
    }
    
    ~PQF(){
        free(ar);
    }
       


    static inline uint64_t bitselect(uint64_t val, int rank) {
        uint64_t i = 1ULL << rank;
        asm("pdep %[val], %[mask], %[val]"
                : [val] "+r" (val)
                : [mask] "r" (i));
        asm("tzcnt %[bit], %[index]"
                : [index] "=r" (i)
                : [bit] "g" (val)
                : "cc");
        return i;
    }


    static inline uint64_t rbitselect(uint64_t val, int rank) {
        uint64_t i = 1ULL << rank;
        asm("pdep %[val], %[mask], %[val]"
                : [val] "+r" (val)
                : [mask] "r" (i));
        asm("lzcnt %[bit], %[index]"
                : [index] "=r" (i)
                : [bit] "g" (val)
                : "cc");
        return i;
    }


    
    inline  int lh2(  int pl, unsigned int filter){
        unsigned long long  int fh =  ( hs2[filter]) ;
        return (int)(fh^pl) ;
    }
    
    inline unsigned  long long int h2(unsigned  long long int pl, unsigned int filter){
        unsigned long long  int fh =  ( hs[filter]) ;
        return (fh^pl) ;
    }
    
    
    inline    unsigned long int hash1 ( unsigned long int ks) { // Bit mix from MurmurHash64/CLHash
        ks ^= ks >> 33;
        ks *= 0xff51afd7ed558ccdULL;
        ks ^= ks >> 33;
        ks *= 0xc4ceb9fe1a85ec53ULL;
        ks ^= ks >> 33;
        return ks;
    }
    
    inline int RS_delete(uint64_t key){
        unsigned long int h = hash1(key);
        unsigned line = h&(num_lines-1);
        h = h>>log_num_lines;
        unsigned place = (h%64);
        h =  h>>6;
        int value = (h&255);
        unsigned char * b =ar + line*64;
        unsigned  long int * is_occ =   (unsigned long int *)(b);
        unsigned  long int * runend =    (unsigned long int *)(b + 8 );
        unsigned  long int extra = 0;
        int last_index = (int)__builtin_popcountll(*is_occ);
        long int last_rank = 63 - rbitselect(rm&*runend , last_index -1 );
        if (last_rank == -1) return 0;
        int index_prev = (int)__builtin_popcountll(*is_occ&((1UL<<place)-1));
        long int rank_prev = 63-rbitselect(rm&*runend , index_prev -1 );
        long int rank_curr = 63-rbitselect(rm&*runend , index_prev);
        int num_pres = (int)(rank_curr - rank_prev);
        int victim = -1;
        for (int i= 0; i < num_pres; i++){
            unsigned char *  src = b + 1 +  16 + rank_prev +  i;
            if ( value == *(src)){victim = i; break;}
        }
        if(victim == -1) return 0;
        int spot_to_del = (int)( rank_prev +  victim);
        unsigned char *  src = b + 1 +  16 + spot_to_del;
        memmove(src, src+1,  last_rank-(1+ spot_to_del) );
        unsigned long int run_mask = (1UL<<(spot_to_del +1 ))-1;
        unsigned long int  top_mask = rm&((~run_mask)>>1);
        if ( num_pres == 1){
            *is_occ &= ~(1UL<<place);  // clear the bit
        }
        *runend =  ( (*runend&~rm) |   (*runend&run_mask) | ( ((rm&*runend)>>1) & top_mask)) & ~extra;
        return 1;
    }

    
    inline    int RS_remove_l(int line, int &place, int &value, int depth){
        //Set Line specific variables
        unsigned char * b =ar + line*64;
        unsigned  long int * is_occ =   (unsigned long int *)(b);
        unsigned  long int * runend =    (unsigned long int *)(b + 8 );
        unsigned  long int extra = 0;
        //Calculate run specific values
        int last_index = (int)__builtin_popcountll(*is_occ);
        long int last_rank = 63 - rbitselect(rm&*runend , last_index -1 );
        int spot = -1;
        int buc = -1;
        buc = spot;
        //If spot is unset then we try to find a random set isOcc bit, if we overflow to the beginning then it resets to the last bucket slot
        if (spot == -1){
            buc = rand()%64;
            int orig_buc = buc;
            while (buc > 0 && !(1 &(*is_occ>>buc))) buc--;
         //   buc_diff += (orig_buc - buc);
            if (buc == 0){
                buc = 63;
                while (buc > 0 && !(1 &(*is_occ>>buc))) buc--;
            //    buc_diff += (64 - buc);
            }
        }
        place = buc;
        // Calculate run indices for the selected run
        int index_prev = (int)__builtin_popcountll(*is_occ&((1UL<<buc)-1));
        long int rank_prev = 63-rbitselect(rm&*runend , index_prev -1 );
        long int rank_curr = 63-rbitselect(rm&*runend , index_prev);
        int num_pres = (int)(rank_curr - rank_prev);
        //Calculate which fingerprint will be taken from this run
        int victim = (int)( rand()%num_pres);
        int spot_to_del = (int)( rank_prev +  victim);
        unsigned char *  src = b + 1 +  16 + spot_to_del;
        value = *(src);
        memmove(src, src+1,  last_rank-(1+ spot_to_del) );
        unsigned long int run_mask = (1UL<<(spot_to_del +1 ))-1;
        unsigned long int  top_mask = rm&((~run_mask)>>1);
        if ( num_pres == 1){
            *is_occ &= ~(1UL<<buc);  // clear the bit
        }
        *runend =  ( (*runend&~rm) |   (*runend&run_mask) | ( ((rm&*runend)>>1) & top_mask)) & ~extra;
        return 0;
    }
    
      
    inline    int RS_setc(uint64_t key){
      unsigned long int h = hash1(key);
      unsigned line = h&(num_lines-1);
      h = h>>log_num_lines;
      unsigned place = (h%64);
      h =  h>>6;
      int value = (h&255);
      return RS_setc_l(line,place, value, 0);
    }

    inline    int RS_setc_l(int line, int place, int value, int depth){
      //Sets preliminary values for the line
        unsigned char * b =ar + line*64;
        unsigned  long int * is_occ =   (unsigned long int *)(b);
        unsigned  long int * runend =    (unsigned long int *)(b + 8 );
        unsigned  long int extra = 0;
        //Calculates run limiters
        int index_prev = (int)__builtin_popcountll(*is_occ&((1UL<<place)-1));
        long int rank_prev = 63-rbitselect(rm&*runend , index_prev -1 );
        int last_index = (int)__builtin_popcountll(*is_occ);
        long int last_rank = 63 - rbitselect(rm&*runend , last_index -1 );
        //If we do not have more space
        if (last_rank > 46){ // we have 48 slots
            if (depth > 256){ cout << "fail on depth " << depth << endl;  return -1; }
            //  cucks++;
            //Remove one item from the line and insert our new value into the filter
            int pl, vl;
            RS_remove_l( line, pl,vl, depth);
            RS_setc_l( line,  place,  value,  depth+1); // This better not fail, as we just deleted something.
            //Calculate the new line and new place for the removed value that we are displacing
            unsigned long int both = (pl<<log_num_lines)| line;
            unsigned long int new_both = h2(both, vl);
            int new_line  = (int)(new_both&(num_lines-1));
            int new_place = (new_both>>log_num_lines)&63;
            //Set the proper OTA bit for the displaced value and insert it in the proper place
            unsigned short * OTA = (unsigned short *)(b +14);
            *OTA |= (1<< (vl%16));
            return RS_setc_l( new_line,  new_place,  vl,  depth+1);
        }
        //If we do have space
        unsigned char *  src = b +  16 + rank_prev +1;
        unsigned long int run_mask = (1UL<<(rank_prev+1))-1;
        unsigned long int  top_mask = rm&((~run_mask)<<1);
        if (!(1&(*is_occ>>place))){
            extra = (1UL << (1+rank_prev));
            *is_occ |= (1UL<<place);
        } else {
            if (*src == value ) return 0;
        }
        memmove(src+1, src, last_rank-rank_prev );
        *src = value;
        *runend =  (*runend&~rm) |  (*runend&run_mask) | ( (*runend<<1) & top_mask) | extra;
        return 0;
    }
    
    inline void add_many( std::vector<uint64_t >& keys,
         std::vector<int>& status,  uint64_t num_keys) {
           for(int i = 0; i < num_keys; i += batch_size){
               std::array<uint64_t, batch_size> hs;
               std::array<uint32_t, batch_size> lines;
               std::array<uint32_t, batch_size> values;
               std::array<uint32_t, batch_size> places;
               std::array<int, batch_size> rank_prevs;
               std::array<bool, batch_size> easy;
               std::array<int, batch_size> last_ranks;
               std::array<unsigned char * , batch_size> bs;
               std::array<uint64_t *, batch_size> is_occs;
               std::array<uint64_t *, batch_size> runends;
               std::array<unsigned short *, batch_size> OTAs;
               std::array<unsigned char * , batch_size> srcs;
                std::array<uint64_t, batch_size> pres;
               for(int j = 0; j < batch_size; j++){
                   hs[j] = hash1(keys[i + j]);
                   pres[j] = 0;
                   easy[j] = false;
               }
               for(int j = 0; j < batch_size; j++){
                   lines[j] = (hs[j]&(num_lines-1));
                   places[j] = ((hs[j]>>log_num_lines)&63);
                   values[j] = (hs[j]>>(6+log_num_lines))&255;
                // }
             //  for(int j = 0; j < batch_size; j++){
                   bs[j] = ar + (lines[j]<<6);
                   is_occs[j] =   (uint64_t *)bs[j];
                   runends[j] =   (uint64_t *)(bs[j] + 8 );
                   OTAs[j] = (unsigned short *)(bs[j] +14);
              // }
              // for(int j = 0; j < batch_size; j++){
                   std::array<int, batch_size> index_prevs;
                    std::array<int, batch_size> last_indexes;
                   index_prevs[j] = (int)__builtin_popcountll((*is_occs[j])&((1UL<<places[j])-1));
                   rank_prevs[j]  = 63-(int)rbitselect((rm&*runends[j]) , index_prevs[j] -1 );
                   last_indexes[j] = (int)__builtin_popcountll((*is_occs[j]));
                   last_ranks[j] = 63 - (int)rbitselect(rm&(*runends[j]) , last_indexes[j] -1 );
               }
               for(int j = 0; j < batch_size; j++){
                   if  (__builtin_expect( ( !(1&( pres[ 127&(lines[j]>>8) ]>>(lines[j]&63)))) && last_ranks[j] < 47 ,1)){
                       easy[j] = true;
                       pres[ 127&(lines[j]>>8)] |= (1UL<<(lines[j]&63));
                       srcs[j] = bs[j] +  16 + rank_prevs[j] +1;
                        memmove(srcs[j]+1, srcs[j], last_ranks[j]-rank_prevs[j] );
                       *srcs[j] = values[j];
                       std::array<uint64_t, batch_size> run_masks;
                       std::array<uint64_t, batch_size> top_masks;
                       std::array<uint64_t, batch_size> extras;
                       run_masks[j] = (1UL<<(rank_prevs[j]+1))-1;
                       top_masks[j] = rm&((~run_masks[j])<<1);
                       extras[j] =0;
                       if (!(1&(*(is_occs[j])>>places[j]))){
                           extras[j] = (1UL << (1+rank_prevs[j]));
                           *(is_occs[j]) |= (1UL<<places[j]);
                       }
                       *runends[j] =  (*runends[j]&~rm) |  (*runends[j]&run_masks[j]) | ( (*runends[j]<<1) & top_masks[j]) | extras[j];
                       status[i+j] = 1;
                   }
               }
               
               for(int j = 0; j < batch_size; j++){
                   if(__builtin_expect(  (!easy[j]),0)){
                       std::array<int, batch_size> index_prevs;
                       std::array<int, batch_size> last_indexes;
                       
                       index_prevs[j] = (int)__builtin_popcountll((*is_occs[j])&((1UL<<places[j])-1));
                       rank_prevs[j]  = 63-(int)rbitselect((rm&*runends[j]) , index_prevs[j] -1 );
                       last_indexes[j] = (int)__builtin_popcountll((*is_occs[j]));
                       last_ranks[j] = 63 - (int)rbitselect(rm&(*runends[j]) , last_indexes[j] -1 );
                       if (__builtin_expect( ( last_ranks[j] > 46),0)){
                           int pl, vl;
                           RS_remove_l( (int)(lines[j]), pl,vl, 0 );
                           int r = RS_setc_l( lines[j],  places[j],  values[j],  1); // This better not fail, as we just deleted something.
                         unsigned long int both = (pl<<log_num_lines)| lines[j];
                           unsigned long int new_both = h2(both, vl);
                           int new_line  = (int)(new_both&(num_lines-1));
                           int new_place = (new_both>>log_num_lines)&63;
                            *(OTAs[j]) |= (1<< (vl%16));
                           r = RS_setc_l( new_line,  new_place,  vl,  1);
                           if (r != 0) cout << "set failed " << endl;
                           status[i+j] = (r == 0);
                       } else{
                           srcs[j] = bs[j] +  16 + rank_prevs[j] +1;
                           memmove(srcs[j]+1, srcs[j], last_ranks[j]-rank_prevs[j] );
                           *srcs[j] = values[j];
                           std::array<uint64_t, batch_size> run_masks;
                           std::array<uint64_t, batch_size> top_masks;
                           std::array<uint64_t, batch_size> extras;
                           run_masks[j] = (1UL<<(rank_prevs[j]+1))-1;
                           top_masks[j] = rm&((~run_masks[j])<<1);
                           extras[j] =0;
                           if (!(1&(*(is_occs[j])>>places[j]))){
                               extras[j] = (1UL << (1+rank_prevs[j]));
                               *(is_occs[j]) |= (1UL<<places[j]);
                           }
                           *runends[j] =  (*runends[j]&~rm) |  (*runends[j]&run_masks[j]) | ( (*runends[j]<<1) & top_masks[j]) | extras[j];
                           status[i+j] = 1;
                       }
                   }
               }
           }
    }
    

   //Top Level Search Function, takes a key and returns 1 if the key is present or 0 otherwise
    inline int RS_getc(uint64_t key){
        unsigned long int h = hash1(key);
        unsigned line = (h&(num_lines-1)); //Using the hash to compute the proper cache line
        h = h>>log_num_lines;
        unsigned place = (h&63); //Using the hash to compute the location within the cache line
        h =  h>>6;
        int value = (h&255);  //Using the hash to compute the fingerprint that we're searching for
        return RS_get_l(line, place, value, false);
    }
    
    //Cache line level search function
    inline int RS_get_l(int line, int place, int value, bool second){
      //Assigning the base variables for the cache line
        unsigned char * b = ar + (line<<6); 
        unsigned long int  is_occ =   *(unsigned long int *)b;
        unsigned long int  runend =    rm & *(unsigned long int *)(b + 8 );
        unsigned short OTA = *(unsigned short *)(b +14);
        //If the proper isOcc bit is unset and there is no bit set for the fingerprint in the OTA then the value is not present
        if (!(1&(is_occ>>place)) &&  !((OTA>>(value%16))&1)){
            return 0;
        }
        //Calculate the beginning and end of the targeted run
        int index_prev = (int)__builtin_popcountll(is_occ&((1UL<<place)-1));
        long int rank_prev  = 63-rbitselect(runend , index_prev -1 );
        int rank_curr = (int)bitselect(runend , index_prev  );
        int num_elements = (int)(rank_curr - rank_prev) ;
        unsigned char * bb = b+ 17 + rank_prev;
         //Check the run for the targeted value
        for ( int j =  0; j <  num_elements; j++){
            if (bb[j] == value ) return 1;
        }
        //If this is the first location and there is an OTA bit for the fingerprint then calculate the values for the other line and recursively call
        if (!second && ((OTA>>(value%16))&1)){
         unsigned long int both = (place<<log_num_lines)| line;
            unsigned long int new_both = h2(both, value);
            int new_line  = (int)(new_both&(num_lines-1));
            int new_place = (new_both>>log_num_lines)&63;
            return RS_get_l(new_line, new_place, value, true);
        }
        return 0;
    }
    
     inline void likely_contains_many( std::vector<uint64_t >& keys,
       std::vector<bool>& status,  uint64_t num_keys) {
         for(int i = 0; i < num_keys; i += batch_size){
             std::array<uint64_t, batch_size> hs;
             std::array<uint64_t, batch_size> lines;
             std::array<uint32_t, batch_size> values;
             std::array<uint64_t, batch_size> places;
             std::array<int, batch_size> rank_prevs;
             std::array<int, batch_size> rank_currs;
             std::array<unsigned char * , batch_size> bs;
             std::array<uint64_t, batch_size> is_occs;
             std::array<uint64_t, batch_size> runends;
             std::array<unsigned short, batch_size> OTAs;
             std::array<bool, batch_size> over;
             std::array<int, batch_size> num_elements;
             std::array< unsigned char *, batch_size> bbs;
             for(int j = 0; j < batch_size; j++){
                 hs[j] = hash1(keys[i + j]);
             }
             for(int j = 0; j < batch_size; j++){
                 lines[j] = (hs[j]&(num_lines-1));
                 places[j] = ((hs[j]>>log_num_lines)&63);
                 values[j] = (hs[j]>>(6+log_num_lines))&255;
             }
             for(int j = 0; j < batch_size; j++){
                 bs[j] = ar + (lines[j]<<6);
                 is_occs[j] =   *(unsigned long int *)bs[j];
                 runends[j] =    rm & *(unsigned long int *)(bs[j] + 8 );
                 OTAs[j] = *(unsigned short *)(bs[j] +14);
                 over[j] = false;
                 status[i+j] = false;
             }
             
             for(int j = 0; j < batch_size; j++){
                 if (!(1&(is_occs[j]>>places[j]))
                     && OTAs[j] == 0
                     ){
                     over[j] = true;
                 }else{
                     std::array<int, batch_size> index_prevs;
                     index_prevs[j] = (int)__builtin_popcountll(is_occs[j]&((1UL<<places[j])-1));
                     rank_prevs[j]  = 63-(int)rbitselect(runends[j] , index_prevs[j] -1 );
                     rank_currs[j] = (int)bitselect(runends[j] , index_prevs[j]  );
                     num_elements[j] = (int)(rank_currs[j] - rank_prevs[j]) ;
                     bbs[j] = bs[j]+ 17 + rank_prevs[j];
                     for ( int k =  0; k <  num_elements[j]; k++){
                         if (bbs[j][k] == values[j] ) {
                             status[i+j] = true;
                             over[j] = true;}
                     }
                 }
             }
             for(int j = 0; j < batch_size; j++){
                 if (!over[j]) {
                      status[i+j] = false;
                     if ( ((OTAs[j]>>(values[j]%16))&1)){
                           unsigned long int both = (places[j]<<log_num_lines)| lines[j];
                         unsigned long int new_both = h2(both, values[j]);
                         int new_line  = (int)(new_both&(num_lines-1));
                         int new_place = (new_both>>log_num_lines)&63;
                          bool found  = (1 == RS_get_l(new_line, new_place, values[j], true));
                         status[i+j] = found;
                         
                     }
                 }
             }
         }
     }

   
    
    
};


#endif /* PQF_h */


/*
 
 
 
 To make these count, we change the following.
 
 When we insert, we check if the value is already present. If it is, we check for a counter and increment, otherwise we add a counter.  If it is not present, we add as normal.  The tricky part is adding a counter.  This moves everything down a byte in the FSA and shifts the run_end as usual.
 
 In lookup we look at the next value, and see if it is a counter, and return the count.
 
 How we encode counters:
 
 A single item is stored as a byte, as usual, but we keep the entire run in order, so we insert in the correct spot, not the first spot as we do now.

 After an item, we keep a count.  If the first byte of the count is smaller tha the item, we recognize it as a count, as fingerpritns should only increase.
 
 If the count is large than the fingerprint, we repeat the fingerprint, a sign that a count will follow.
 
 Counts are encoded like the are in Ethereum.  Value between 0 and 250 as stored as bytes.  Values higher than that mean 251 - the next two bytes are a value, 252: next three bytes are a value, etc.
 
 The cases are:
 1.  Item not present, insert as usual, but get order correct.
 2. Item present, with no count.  Two subcases, item is 0 (or 1), in which case we repeat the fingerprint, and follow with by the byte 1.  This inserts two bytes.
    If the item is >1, we insert 1 byte, the count 1, which follows the item.
 3. The item is present with a count < 250.  We increment the count.
 4. The item is present with a count = 250.  We change 250 to 251, and add 2 bytes,
 5. The item is present with a count > 250.  We increment the count, and if we pass 65k, as change 251 to 252, and add a byte.
 
 
 If we need to evict a value, we choose one that does not have a count.
 
 */
