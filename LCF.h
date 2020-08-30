//
//  LCF.h
//  LCF
//


#ifndef LCF_h
#define LCF_h


//static  //256;


class LCF {
   static  constexpr int batch_size = 128;
    
    static  constexpr int cacheline_size = 64;
    static constexpr int fingerprint_max = 256;
     static constexpr int fingerprint_mask  = 255;
    static  constexpr int last_place_in_cache_line = 63;
    
    static  constexpr int number_of_virtual_slots = 96;
   
    static constexpr int start_of_fingerprint_storage = 14;
        static constexpr int start_of_OTA = 12;
    
    static constexpr int offset_shift = 6;
    static constexpr int offset_mask = 63;
    
    static constexpr uint64_t mask_for_slot_array = 0xffffffff;
    
    static constexpr int tot_entries_thres = 49;
    
     static constexpr int  OTA_size = 16;
    
    static constexpr int place_mask = 255;
    static constexpr int place_shift = 10;
       static constexpr int line_shift = 17;
    
    
public:
    
    unsigned char * ar1; //
    unsigned char * ar_end; //
 
    unsigned long long int  hs2[fingerprint_max];
    unsigned long long int  hs3[fingerprint_max];
    long  num_lines;
    long num_lines1;
    int lg_num_lines;
    
    LCF( long size){
        
        long sz = size; // one cache line
        if (posix_memalign((void **)&ar1,cacheline_size, sz * cacheline_size ) != 0) {
            cout << "error on memalign" << endl;
            exit(1);
        }
        memset (ar1,1,sz*cacheline_size);  // fault the memory in?
        bzero(ar1,sz*cacheline_size);
        num_lines = sz;
        num_lines1 = num_lines-1;
        lg_num_lines = (int)log2(num_lines);
       
        // hs2 is for lines
        for (int i =0;i < fingerprint_max; i++){
            hs2[i] = ( (unsigned  int)hash1( rand() +1 )&(0x1fff  ));
            while ((hs2[i]&(0x1fff )) == 0 )
                hs2[i] = ( (unsigned  int)hash1( rand() +1 )&(  0x1fff))   ;
        }
        //hs3 is for places (0--95)
        for (int i =0;i < fingerprint_max; i++){
            hs3[i] = ( (unsigned  int)hash1(rand() )&( 0x3f  ));
            while ((hs3[i]&(0x3f))  == 0  )
                hs3[i] = ( (unsigned  int)hash1( rand() +1 )&( 0x3f))   ;
        }
       // print = false;
        
        
    }
    
    ~LCF(){
        free(ar1);
    }
    
    
    inline     int PCF_get1(unsigned char * ar, int n){
        return ar[start_of_fingerprint_storage + n];}
    
    
    inline    int PCF_set1(unsigned char * ar, int n, int v){
        return ar[start_of_fingerprint_storage + n] = v;
    }
    
    inline int compute_entries_to_left(unsigned long long * fca, int step,  unsigned long m0){
         return (int)__builtin_popcountll( (fca[step])&m0) + step * (int)__builtin_popcountll(*fca);
     }
    
    inline unsigned int CF1_get(int bl,int buc){
        unsigned char * b = ar1 + (bl * cacheline_size);
        unsigned long long * fca = (unsigned long long *)b ;
        unsigned long m0, bit_to_set;
        int off,step;
        compute_offsets(m0,  bit_to_set,off, step, b, buc );
        if (! ((fca[step]) &  bit_to_set)){
            return -1;
        }
        int entries_to_left =  compute_entries_to_left(fca,  step,  m0);
        return  PCF_get1(b,entries_to_left);
    }
    
    
    inline void compute_offsets( unsigned long &m0,  unsigned long &bit_to_set,int &off, int &step, unsigned char * b, int place ){
        step = place>>offset_shift;
        off = place&offset_mask;
        bit_to_set = (1UL<<off);
        m0 = ( bit_to_set -1);
    }
    
    inline int compute_tot_entries(unsigned long long * fca){
        int one_bits_set = (int)__builtin_popcountll(*fca );
        int two_bits_set = (int)__builtin_popcountll(fca[1] & mask_for_slot_array);
        return one_bits_set + two_bits_set ;
    }
    
    inline unsigned int CF1_set(int bl,int buc, int v){
        unsigned char * b = ar1 + (bl * cacheline_size);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        unsigned long m0, bit_to_set;
        int off,step;
        compute_offsets(m0,  bit_to_set,off, step, b, buc );
        int entries_to_left =  compute_entries_to_left(fca,  step,  m0);
        int tot_entries = compute_tot_entries(fca);
        if (tot_entries > tot_entries_thres ){
            return -1;
        }
        if ( !((fca[step]) &  bit_to_set) ){
            memmove( b + start_of_fingerprint_storage + 1 +  entries_to_left,  b + start_of_fingerprint_storage + entries_to_left, tot_entries-entries_to_left);
        }
        PCF_set1(b,  entries_to_left  , v);
        fca[step]  |=    bit_to_set ;
        return 0;
    }
    
    inline unsigned int CF2_set(int bl,int buc, int v){ // return -1 for full, 1 for success, 0 for already occupied
        unsigned char * b = ar1 + (bl * cacheline_size);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        unsigned long m0, bit_to_set;
        int off,step;
        compute_offsets(m0,  bit_to_set,off, step, b, buc );
        if ( ((fca[step]) &  bit_to_set) ){
            return 0;
        }
        int entries_to_left =  compute_entries_to_left(fca,  step,  m0);
        int tot_entries =  compute_tot_entries(fca);
        if (tot_entries > tot_entries_thres ){
            return -1;
        }
        memmove( b + start_of_fingerprint_storage + 1 +  entries_to_left,  b + start_of_fingerprint_storage + entries_to_left, tot_entries-entries_to_left);
        PCF_set1(b,  entries_to_left  , v);
        fca[step]  |=    bit_to_set ;
        return 1;
    }
     
    
    inline    unsigned long int hash1 (unsigned long int ks) { // Bit mix from MurmurHash64/CLHash
        ks ^= ks >> 33;
        ks *= 0xff51afd7ed558ccdULL;
        ks ^= ks >> 33;
        ks *= 0xc4ceb9fe1a85ec53ULL;
        ks ^= ks >> 33;
        return ks;
    }
 
    // This maps local places to local places.
    inline  int h2(  int pl, unsigned int filter){
        long one = 1;
        long long  int   offset = ((hs3[filter] & 0x3f) + 1) | one; // The offset is at most 63.
        offset = (pl & one) ? offset : -offset;
        int64_t output = pl + offset;
        if(output < 0){output += number_of_virtual_slots;   }
        if(output >= number_of_virtual_slots){   output -= number_of_virtual_slots; }
        return (int)output;
    }
    
    
    // Warning, if 0x1fff is greater than num_lines, bad things can happen. 8192
    // makes a new line from an old, arbitrary num_lines.
    inline  int lh2(  int pl, unsigned int filter){
        long one = 1;
        long long  int   offset = ((hs2[filter] & 0x1fff) + (1)) | one;
        offset = (pl & one) ? offset : -offset;
        int64_t output = static_cast<int64_t>(pl) + offset;
        if(output < 0){            output += num_lines;  }
        if(static_cast<uint64_t>(output) >= num_lines){     output -= num_lines;}
        return (int)output;
    }
    

    // makes a new place 0..95 from an old.  This maps  a place to a place in the new block.
    inline  int lh3(  int pl, unsigned int filter){
        return pl;
        }
   
    
    inline  int agg_hget_l_alt( int bl, int place, int place1, int filter){
        unsigned char * b = ar1 + (bl * cacheline_size);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        unsigned long m0, bit_to_set;
        int step,off;
        compute_offsets(m0,  bit_to_set,off, step, b, place );
        if (!((((unsigned long long *)b)[step]) &  bit_to_set) ){return 0; }
        int entries_to_left = compute_entries_to_left(fca,  step,  m0);
        if ( __builtin_expect( (filter ==  PCF_get1(b,entries_to_left)) ,0))  return 1;
        compute_offsets(m0,  bit_to_set,off, step, b, place1 );
        if (!((((unsigned long long *)b)[step]) &  bit_to_set)){  return 0; }
        int  entries_to_left1 = compute_entries_to_left(fca,  step,  m0);
        if ( __builtin_expect( (filter ==  PCF_get1(b,entries_to_left1)),0)) return 1;
        return 0;
        
    }
      
    inline void compute_new_block(int place, int place1,int filter, int bl, int &new_bl, int &new_place, int &new_place1){
        int lower = (place<place1)?place:place1;
        new_bl = lh2(bl,filter);
        new_place = lh3(lower,filter);
        new_place1 = h2(new_place, filter);   }
    
    inline  int agg_hget(uint64_t key){
        unsigned long  int h = (unsigned long int) hash1(key);
        int filter = ( h & fingerprint_mask);
        int place =  (h>>place_shift)&place_mask;
        place = (place*number_of_virtual_slots)>>8;
        int bl = (int)((h>>line_shift)&(num_lines1));
        unsigned char * b = ar1 + (bl * cacheline_size);
              unsigned long long *  fca = ( unsigned long long  *)b ;
        unsigned long m0, bit_to_set;
        int off,step;
        compute_offsets(m0,  bit_to_set,off, step, b, place );
        unsigned short OTA = *(unsigned short *)(b + start_of_OTA);
        if (OTA == 0){
            if (!((((unsigned long long *)b)[step]) &  bit_to_set) ){return 0; }
            if (filter ==  PCF_get1(b,compute_entries_to_left(fca,  step,  m0))) return 1;
            compute_offsets(m0,  bit_to_set,off, step, b,  h2(place,filter));
            if (! ((((unsigned long long *)b)[step]) &  bit_to_set)){  return 0; }
            if (filter ==  PCF_get1(b,compute_entries_to_left(fca,  step,  m0))) return 1;
            return 0;
        }
        if (  ((((unsigned long long *)b)[step]) &  bit_to_set) ){
            if (filter ==  PCF_get1(b,compute_entries_to_left(fca,  step,  m0))) return 1;
        }
        int place1 = h2(place,filter);
        compute_offsets(m0,  bit_to_set,off, step, b, place1 );
        if ( ((((unsigned long long *)b)[step]) &  bit_to_set)  ){
            if (filter ==  PCF_get1(b,compute_entries_to_left(fca,  step,  m0))) return 1;
        }
        if  (1&(OTA>>(filter%OTA_size) )){
            int new_bl, new_place, new_place1;
            compute_new_block( place,  place1, filter,  bl, new_bl,new_place, new_place1);
            return hget_l( new_bl,  new_place, new_place1,  filter ,true);
        }
        return 0;
    }
    
    inline void generate_block_and_place( unsigned long  int h, int &filter, int &place,  int &bl ){
        filter = ( h & fingerprint_mask);
        h =h>>place_shift;
        place = h&place_mask;
        place = (place*number_of_virtual_slots)>>8;
        h=(h>>(line_shift-place_shift ));
        bl = (int)(h&(num_lines -1 ));
    }
    
    inline  int hget(uint64_t key){
        int filter, place, bl;
        generate_block_and_place(   hash1(key), filter,  place,   bl );
        int place1 = h2(place,filter);
        return hget_l(bl, place, place1, filter, false);
    }
    
    inline  int hget_l( int bl, int place, int place1, int filter , bool secondary){
        int r = CF1_get(bl,place);
        if (r == filter){
            return 1;
        }
        unsigned short OTA = *(unsigned short *)(ar1 + (bl * cacheline_size) + start_of_OTA);
        if ( r == -1 &&   !secondary && OTA == 0 ){
            return 0;
        }
        if (( CF1_get(bl,place1 )) == filter){
            return 1;
        }
        if  ( !secondary && 1&(OTA>>(filter%OTA_size) )){
            int new_bl, new_place, new_place1;
            compute_new_block( place,  place1, filter,  bl, new_bl,new_place, new_place1);
            return hget_l( new_bl,  new_place, new_place1,  filter ,true);
        }
        return 0;
    }
    
    
    inline  int  is_room(int bl){
        unsigned char * b = ar1 + (bl * cacheline_size);
        unsigned long long  * fca = ( unsigned long long *)b ;
        int tot_entries = compute_tot_entries(fca);
        return  51 - tot_entries ;
    }
    
    inline bool quick_move(int bl, int place, int filt,  vector< pair < int, int > > &to_move_list ){
        int tries = 0;
        to_move_list.push_back(make_pair(place, filt));
        int pos = place;
        int filter = CF1_get(bl,pos );
        if (filter == -1){ return true;     }
        while(tries++ < 16)
        {
            //moves++;
            int next = h2(pos, filter );
            int new_filter =  CF1_get(bl,next ) ;
            if (new_filter == filter) {  return true; }
            if ( next == place){
                for (auto k: to_move_list){
                    if (place == k.first && filter == k.second){ return true;  }
                }
                to_move_list.push_back(make_pair(next, filter));
                return false;
            }
            to_move_list.push_back(make_pair(next, filter));
            if(new_filter == -1){  // This needs one more space, so may fail.
                if (is_room(bl) > 0 ) { return true;}
                return false;
            }
            filter = new_filter;
            pos = next;
        }
        return false;  // the fall out case.
    }
    
    inline   int cf_move(int bl, int place, int filt, bool early, vector < pair < int ,int > > &v_out){
        //cucks++;
        vector< pair < int, int > > to_move_list;
        bool rdone = quick_move( bl,  place,  filt, to_move_list );
        if (!rdone) {
            to_move_list.clear();
            rdone = quick_move( bl,  h2(place,filt),  filt, to_move_list );
        }
        if ( rdone){
            for (auto y : to_move_list){
                int r =   CF1_set(bl,y.first,y.second);
                if (r == -1){
                    //  We check in quick move that this can never happen.
                    return -1;
                }
            }
            return 1;
        }
        v_out =  to_move_list;
   //     cuckoo_fail++;
        return -2;
    }
    
    
    
    inline  int hput(uint64_t key){
        int filter, place, bl;
        generate_block_and_place(   hash1(key), filter,  place,   bl );
        int place1 = h2(place,filter);
        return  hput_l(bl,place,place1, filter, key, false,0);
        
    }
    
    inline int delete_entry( int bl,int place, int & filter){
        //  removed++;
        unsigned char * b = ar1 + (bl * cacheline_size);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        unsigned long m0, bit_to_set;
        int step, off;
        compute_offsets(m0,  bit_to_set,off, step, b, place );
        int entries_to_left = compute_entries_to_left(fca,  step,  m0);
        int tot_entries = compute_tot_entries(fca);
        filter = PCF_get1(b,  entries_to_left );
        memmove( b + start_of_fingerprint_storage +  entries_to_left,  b + 1+ start_of_fingerprint_storage + entries_to_left, tot_entries-(1+ entries_to_left));
        fca[step]  &=    ~bit_to_set ;
     //   lcucks++;
        return 0;
    }
    
    
    inline int is_present(unsigned char * b, int off){
        unsigned long long *  fca = ( unsigned long long  *)b ;
        unsigned long bit_to_set;
        int step = off>>offset_shift;
        off = off&offset_mask;
        bit_to_set = (1UL<<off);
        return  (0!= (fca[step] &  bit_to_set));
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
    
    
    
    inline      int line_full( int bl,int place, int place1, int filter, uint64_t key, bool early, int depth){
        if( depth > 450){
            cout << "Depth too high " << depth << " bl " << bl << endl;
            return -1;
        }
        //lines_full++;
        unsigned char * b = ar1 + (bl * cacheline_size);
        unsigned short * OTA = (unsigned short *)(ar1 + (bl * cacheline_size) + start_of_OTA);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        int spot = -1;
        int spot1 = -1;
        int first_empty = -1;
        bool emp = false;
        if ( OTA!= 0  && depth < 16
            ){
            //  line_full_ota++;
            for (int k = 0; k < tot_entries_thres; k++){
                if (  ((short)(*OTA>>(b[start_of_fingerprint_storage+k]%OTA_size))&1) != 0 ){
                    int new_bl = lh2(bl, b[start_of_fingerprint_storage+k] );
                    if (ar1[new_bl*cacheline_size + last_place_in_cache_line ] == 0){
                        emp = true;
                        //  found_empty++;
                        int n0 = (int)__builtin_popcountll(fca[0]);
                        spot = ( k >= n0)?(64 + (int)bitselect(fca[1], k-n0)):(int)bitselect(fca[0], k);
                        break;
                    }
                    
                } else {
                    if ( first_empty == -1 && k > depth/16 ){
                        int new_bl = lh2(bl, b[start_of_fingerprint_storage+k] );
                        if (ar1[new_bl*   cacheline_size + last_place_in_cache_line] == 0){
                            first_empty = k;
                            int n0 = (int)__builtin_popcountll(fca[0]);
                            spot1 = ( k >= n0)?(64 + (int)bitselect(fca[1], k-n0)):(int)bitselect(fca[0], k);
                        }
                    }
                }
            }
        }
        //   if (!emp && spot != -1) found_nonempty++;
        if (spot == -1 ){
            //   if ( OTA!= 0 && depth < 16  ) line_full_zero_ota++;
             spot = spot1;
            }
         if (spot == -1 ){
            //   line_full_rand++;
             int buc = rand()%number_of_virtual_slots;
             while(buc > 0 &&  !is_present(b,buc) ){ buc--;  }
             if (buc == 0){
                 buc = number_of_virtual_slots-1;
                 while(buc > 0 &&  !is_present(b,buc) ){ buc --; }
             }
             spot = buc;
             
         }
           int pl = spot;
       int fl1;
     delete_entry(bl,pl, fl1);
        int r =  agg_hput_l( bl, place,  place1,  filter,  key,  false, depth+1);
        if (r == -1) {
            *OTA |= (1<<(fl1%OTA_size));
            return -1;
        }
        int pl1 = h2(pl,fl1);
        r =   hput_l_alt( bl, pl,  pl1,  fl1,  0,  false, depth+1);
        if (r == -1) {
            *OTA |= (1<<(fl1%OTA_size));
            return -1;
        }
        *OTA |= (1<<(fl1%OTA_size));
      return 1;
        
    }
    
     inline int loop_failure( int bl,int place, int place1, int filter, uint64_t key, bool early, vector < pair< int, int> > v, int depth){
        //loop_fail++;
        int item_no = rand()%(v.size());
        int pl = v[item_no].first;
        int fl1;
        delete_entry(bl,pl, fl1);
        int r =  agg_hput_l( bl, place,  place1,  filter,  key,  false, depth+1);
        unsigned short * OTA = (unsigned short *)(ar1 + (bl * cacheline_size) + start_of_OTA);
        if (r == -1) {
            *OTA |= (1<<(fl1%OTA_size));
            return -1;
        }
        int pl1 =  h2(pl,fl1);
        r = hput_l_alt( bl, pl, pl1,  fl1,  0,  false, depth+1);
        if (r == -1){
            *OTA |= (1<<(fl1%OTA_size));
            return -1;
        }
         *OTA |= (1<<(fl1%OTA_size)); //GGG
        return 1;
        
     }
    
    inline  int hput_l_alt(int bl, int place, int place1, int filter, int key, bool early, int depth){
        int new_bl, new_place, new_place1;
        compute_new_block( place,  place1, filter,  bl, new_bl,new_place, new_place1);
        return hput_l( new_bl,  new_place, new_place1,  filter ,key,true, depth+1);  // why true?  Because this is is a secondary insert?  Maybe this is not used.
        
    }
    
  //  #define ENABLE_PREFETCH
    
   
 
    inline void likely_contains_many( std::vector<uint64_t >& keys,
                                     std::vector<bool>& status,  uint64_t num_keys) {
        for(int i = 0; i < num_keys; i += batch_size){
            std::array<uint64_t, batch_size> hs;
            std::array<uint64_t, batch_size> bls;
            std::array<uint32_t, batch_size> values;
            std::array<int, batch_size> places;
            std::array<int, batch_size> place1s;
            std::array<unsigned char *, batch_size> bs;
            std::array<int, batch_size> offs;
            std::array<uint64_t, batch_size> m0s;
            std::array<uint64_t, batch_size> bit_to_sets;
            std::array<int, batch_size> steps;
            std::array<int, batch_size> entries_to_lefts;
            std::array<bool, batch_size> over;
            std::array<short, batch_size> OTAs;
            for(int j = 0; j < batch_size; j++){
                hs[j] = hash1(keys[i + j]);
            }
            for(int j = 0; j < batch_size; j++){
                values[j] = (hs[j])&fingerprint_mask;
                bls[j] = ((hs[j]>>line_shift)&(num_lines1));
                places[j] = ((hs[j]>>place_shift)&place_mask);
                   places[j] = (places[j]*number_of_virtual_slots)>>8;
             //   places[j] = (places[j]>95)?places[j]-64:places[j];
                }
            for(int j = 0; j < batch_size; j++){
                //offs[j] = (int)places[j];
                bs[j] = ar1 + (bls[j] * cacheline_size);
                _mm_prefetch((char*)( bs[j] ),_MM_HINT_T0);
            }
            for(int j = 0; j < batch_size; j++){
                steps[j] = places[j]>>offset_shift;
                offs[j] = places[j]&offset_mask;
                bit_to_sets[j] = (1UL<<offs[j]);
            }
            for(int j = 0; j < batch_size; j++){
                OTAs[j] = *(unsigned short *)(bs[j] + start_of_OTA);
                over[j] = false;
            }
            for(int j = 0; j < batch_size; j++){
                if ( OTAs[j] == 0){
                    if (!((((unsigned long long *)bs[j])[steps[j]]) &  bit_to_sets[j]) ){
                        //    first_fast_fail++;
                        status[i+j] = false;   over[j] = true; }
                    else{
                        m0s[j] = ( bit_to_sets[j] -1);
                        entries_to_lefts[j] = (int)__builtin_popcountll( (((unsigned long long *)bs[j])[steps[j]])&m0s[j]) +
                        steps[j] *  (int)__builtin_popcountll(((unsigned long long *)bs[j])[0]);
                        if (values[j] ==  PCF_get1(bs[j],entries_to_lefts[j])){
                            //   first_look_success++;
                            status[i+j] = true;   over[j] = true; }
                        else {
                            //    first_look_fail++;
                            place1s[j] = h2(places[j],values[j]);
                            offs[j] = place1s[j];
                            steps[j] = offs[j]>>offset_shift;
                            offs[j] = offs[j]&offset_mask;
                            bit_to_sets[j] = (1UL<<offs[j]);
                            if (! ((((unsigned long long *)bs[j])[steps[j]]) &  bit_to_sets[j])){
                                //  second_fast_fail++;
                                status[i+j] = false;   over[j] = true; }
                            else {
                                m0s[j] = ( bit_to_sets[j] -1);
                                entries_to_lefts[j] = (int)__builtin_popcountll( (((unsigned long long *)bs[j])[steps[j]])&m0s[j]) +
                                steps[j] *   (int)__builtin_popcountll(((unsigned long long *)bs[j])[0]);
                                status[i+j] = (values[j] ==  PCF_get1(bs[j],entries_to_lefts[j]));
                                if (values[j] ==  PCF_get1(bs[j],entries_to_lefts[j])){
                                    // second_look_success++;
                                    
                                }
                                else {
                                    //second_look_fail++;}
                                    over[j] = true;
                                }
                            }
                        }
                    }
                }
            }
            
            for(int j = 0; j < batch_size; j++){
                if (!over[j]){
                    place1s[j]  = h2(places[j] ,values[j]);
                    //  int lower = (places[j] <place1s[j] )?places[j] :place1s[j] ;
                    if (  ((((unsigned long long *)bs[j] )[steps[j] ]) &  bit_to_sets[j] ) ){
                        unsigned long m0s[j] ;
                        m0s[j]  = ( bit_to_sets[j]  -1);
                        entries_to_lefts[j]  = (int)__builtin_popcountll( (((unsigned long long *)bs[j] )[steps[j] ])&m0s[j] ) + steps[j]  * (int)__builtin_popcountll(((unsigned long long *)bs[j] )[0]) ;
                        if (values[j] ==  PCF_get1(bs[j] ,entries_to_lefts[j] )) { status[i+j] = true;   over[j] = true; }
                    }
                    if (!over[j]){
                        offs[j]  = place1s[j] ;
                        steps[j]  = offs[j] >>offset_shift;
                        offs[j]  = offs[j] &offset_mask;
                        bit_to_sets[j]  = (1UL<<offs[j] );
                        if ( ((((unsigned long long *)bs[j] )[steps[j] ]) &  bit_to_sets[j] )  ){
                            m0s[j]  = ( bit_to_sets[j]  -1);
                            entries_to_lefts[j]  = (int)__builtin_popcountll( (((unsigned long long *)bs[j] )[steps[j] ])&m0s[j] ) + steps[j]  * (int)__builtin_popcountll(((unsigned long long *)bs[j] )[0]) ;
                            if (values[j] ==  PCF_get1(bs[j] ,entries_to_lefts[j] )) {  status[i+j] = true;   over[j] = true; }
                        }
                    }
                    if  (!over[j] && (1&(OTAs[j]>>(values[j]%OTA_size) ))){ // GGG
                         int lower = (places[j] <place1s[j] )?places[j]:place1s[j] ;
                        int new_bl = lh2((int)bls[j],values[j]);
                        int new_place = lh3(lower,values[j]);
                          int new_place1 = h2(new_place, values[j]);
                         status[i+j] = hget_l( new_bl,  new_place, new_place1,  values[j] ,true);
                    } else if  (!over[j] ){
                       // ota_miss++;
                        status[i+j] = false;
                    }
                }
            }
            }
    }
    
    inline  int hdelete(uint64_t key){
        int filter, place, bl;
        generate_block_and_place(   hash1(key), filter,  place,   bl );
        int place1 = h2(place,filter);
        return hdelete_l( bl, place,  place1,  filter,  false);
    }
    
    
    inline    int hdelete_l(int bl,int place, int place1, int filter,  bool sec){
        unsigned char * b = ar1 + (bl * cacheline_size);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        unsigned long m0, bit_to_set;
        int off,step;
        compute_offsets(m0,  bit_to_set,off, step, b, place );
        if (__builtin_expect( ((fca[step]) &  bit_to_set),0 )){  // already occupied, try alt position optimize for space available.
            int entries_to_left =  compute_entries_to_left(fca,  step,  m0);
            int tot_entries =  compute_tot_entries(fca);
            if (b[start_of_fingerprint_storage+entries_to_left] == filter) {
                memmove( b + start_of_fingerprint_storage +  entries_to_left,  b + start_of_fingerprint_storage + 1+ entries_to_left, tot_entries-entries_to_left);
                fca[step]  &=    ~bit_to_set ;
                return 1;
            }
        }
        compute_offsets(m0,  bit_to_set,off, step, b, place1 );
        if ( __builtin_expect(  ((fca[step]) &  bit_to_set),0 )){ // optimize for space available.
            int entries_to_left =  compute_entries_to_left(fca,  step,  m0);
            int tot_entries =  compute_tot_entries(fca);
            if (b[start_of_fingerprint_storage+entries_to_left] == filter){
                memmove( b + start_of_fingerprint_storage +  entries_to_left,  b + start_of_fingerprint_storage + 1+ entries_to_left, tot_entries-entries_to_left);
                fca[step]  &=  ~bit_to_set ;
                return 1;
            }
        }
        if(sec) return 0;
        int new_bl, new_place, new_place1;
        compute_new_block( place,  place1, filter,  bl, new_bl,new_place, new_place1);
        return hdelete_l( new_bl, new_place,  new_place1,  filter,  true);
        
    }
    
    
    inline  int agg_hput(uint64_t key){
        int filter, place, bl;
        generate_block_and_place(   hash1(key), filter,  place,   bl );
        int place1 = h2(place,filter);
        return  agg_hput_l(bl,place,place1, filter, key, false,0);
        
    }
    
    
    inline    int agg_hput_l(int bl,int place, int place1, int filter, uint64_t key, bool early, int depth){
        unsigned char * b = ar1 + (bl * cacheline_size);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        unsigned long m0, bit_to_set;
        int off,step;
        compute_offsets(m0,  bit_to_set,off, step, b, place );
        if (__builtin_expect( ((fca[step]) &  bit_to_set),0 )){  // already occupied, try alt position optimize for space available.
            int entries_to_left =  compute_entries_to_left(fca,  step,  m0);
            if (b[start_of_fingerprint_storage+entries_to_left] == filter) return 1;
              compute_offsets(m0,  bit_to_set,off, step, b, place1 );
            if ( __builtin_expect(  ((fca[step]) &  bit_to_set),0 )){ // optimize for space available.
                int entries_to_left =  compute_entries_to_left(fca,  step,  m0);
                if (b[start_of_fingerprint_storage+entries_to_left] == filter) return 1;
                int tot_entries =  compute_tot_entries(fca);
                if (__builtin_expect(  (tot_entries > tot_entries_thres),0)){ return line_full(  bl, place,  place1,  filter,  key,  false, depth); } // optimize presuming space
                vector < pair < int ,int > > v;
                int r1 = cf_move(bl,place, filter, early,v);
                if (r1 == 1) return 1;
                if (r1 ==-1){ return line_full(  bl, place,  place1,  filter,  key,  false, depth);    }
                if (r1 ==-2){   return   loop_failure(  bl, place,  place1,  filter,  key,  false,  v, depth);   }
                return -1;
            }
        }
        int entries_to_left =  compute_entries_to_left(fca,  step,  m0);
        int tot_entries =  compute_tot_entries(fca);
        if (__builtin_expect( (tot_entries > tot_entries_thres ),0)){  // optimize presuming space
            return line_full(  bl, place,  place1,  filter,  key,  false, depth);
        } else {
            memmove( b + start_of_fingerprint_storage + 1 +  entries_to_left,  b + start_of_fingerprint_storage + entries_to_left, tot_entries-entries_to_left);
            PCF_set1(b,  entries_to_left  , filter);
            fca[step]  |=    bit_to_set ;
            return 1;
        }
    }
    
    
    inline void add_many( std::vector<uint64_t >& keys,
                         std::vector<bool>& status,  uint64_t num_keys) {
        for(int i = 0; i < num_keys; i += batch_size){
            std::array<uint64_t, batch_size> hs;
            std::array<int, batch_size> bls;
            std::array<uint32_t, batch_size> values;
            std::array<int, batch_size> places;
            std::array<int, batch_size> place1s;
            std::array<unsigned char *, batch_size> bs;
            std::array<uint64_t *, batch_size> fcas;
            std::array<int, batch_size> offs;
            std::array<uint64_t, batch_size> m0s;
            std::array<uint64_t, batch_size> bit_to_sets;
            std::array<int, batch_size> steps;
            std::array<int, batch_size> entries_to_lefts;
            std::array<int, batch_size> tot_entries;
            std::array<bool, batch_size> over;
            //   std::array<short, batch_size> OTAs;
            std::array<uint64_t, batch_size> pres;
            std::array<bool, batch_size> easy;
            //    std::array<bool, batch_size> first_test;
            for(int j = 0; j < batch_size; j++){
                hs[j] = hash1(keys[i + j]);
                pres[j] = 0;
                easy[j] = false;
            }
            for(int j = 0; j < batch_size; j++){
                values[j] = (hs[j])&fingerprint_mask;
                bls[j] = (int)((hs[j]>>line_shift)&(num_lines1));
                places[j] = ((hs[j]>>place_shift)&place_mask);
                   places[j] = (places[j]*number_of_virtual_slots)>>8;
//                places[j] = (places[j]>95)?places[j]-64:places[j];
                             
                offs[j] = places[j];
                bs[j] = ar1 + (bls[j] * cacheline_size);
                fcas[j] = ( uint64_t *)bs[j] ;
                steps[j] = (offs[j]/64);
                offs[j] = offs[j]&offset_mask;
                bit_to_sets[j] = (1UL<<offs[j]);
                m0s[j] = ( bit_to_sets[j] -1);
                over[j] = false;
                entries_to_lefts[j]  = (int)__builtin_popcountll( (fcas[j][steps[j] ])&m0s[j] ) + steps[j]  * (int)__builtin_popcountll(fcas[j][0]) ;
                tot_entries[j] = (int)__builtin_popcountll( fcas[j][0]) + (int)__builtin_popcountll(fcas[j][1] & mask_for_slot_array ) ;
                
            }
            for(int j = 0; j < batch_size; j++){
                if  ( __builtin_expect( ( ( !(1&( pres[ 127&(bls[j]>>8) ]>>(bls[j]&offset_mask))))
                                         && !  ((fcas[j][steps[j]]) &  bit_to_sets[j]) && tot_entries[j] <= 45 )   ,1)){
                    easy[j] = true;
                    // HHH no check that it is already present!!!!
                    // cout << "before  " <<  std::hex << pres[ 127&(places[j]>>8)] << std::dec << endl;
                    pres[ 127&(bls[j]>>8)] |= (1UL<<(bls[j]&offset_mask));
                    if (__builtin_expect( (bs[j][start_of_fingerprint_storage+entries_to_lefts[j]] == values[j]),0)) { status[i+j] =  true ;continue; }
                    memmove( bs[j] + start_of_fingerprint_storage + 1 +  entries_to_lefts[j] ,  bs[j] + start_of_fingerprint_storage + entries_to_lefts[j] , tot_entries[j] -entries_to_lefts[j] );
                    PCF_set1(bs[j],  entries_to_lefts[j]   , values[j]);
                    fcas[j] [steps[j] ]  |=    bit_to_sets[j]  ;
                    status[i+j] = true;
                }
            }
            for(int j = 0; j < batch_size; j++){
                if (__builtin_expect( (!easy[j]),0)){
                    place1s[j] = h2(places[j],values[j]);
                    status[i+j] = (1 == hput_l(bls[j], places[j],place1s[j], values[j],0,false, 0));
                    continue;
                    if (__builtin_expect( ((fcas[j][steps[j]]) &  bit_to_sets[j]),0 )){  // already occupied, try alt position optimize for space available.
                        //TJC
                        tot_entries[j] = (int)__builtin_popcountll( fcas[j][0]) + (int)__builtin_popcountll(fcas[j][1] & 0xfffffff ) ;
                        entries_to_lefts[j] = (int)__builtin_popcountll( (fcas[j][steps[j]])&m0s[j]) + steps[j] * (int)__builtin_popcountll(fcas[j][0]) ;
                        if (bs[j][start_of_fingerprint_storage+entries_to_lefts[j]] == values[j]) { status[i+j] =  true;continue; }
                        //TJC
                        place1s[j] = h2(places[j],values[j]);
                        offs[j] = place1s[j];
                        steps[j] = (offs[j]/64);
                        offs[j] = offs[j]&offset_mask;
                        bit_to_sets[j] = (1UL<<offs[j]);
                        m0s[j] = ( bit_to_sets[j] -1);
                        if ( __builtin_expect(  ((fcas[j][steps[j]]) &  bit_to_sets[j]),0 )){ // optimize for space available.
                            entries_to_lefts[j] = (int)__builtin_popcountll( (fcas[j][steps[j]])&m0s[j]) + steps[j] * (int)__builtin_popcountll(fcas[j][0]) ;
                            if (bs[j][start_of_fingerprint_storage+entries_to_lefts[j]] == values[j])  { status[i+j] =  true; continue; }
                            //TJC
                            tot_entries[j] = (int)__builtin_popcountll( fcas[j][0]) + (int)__builtin_popcountll(fcas[j][1] & mask_for_slot_array) ;
                            if (__builtin_expect(  (tot_entries[j] > tot_entries_thres),0)){ status[i+j] =  1==line_full(  bls[j], places[j],  place1s[j],  values[j],  0,  false, 0); } // optimize presuming space
                            vector < pair < int ,int > > v;
                            int r1 = cf_move(bls[j],places[j], values[j], false,v);
                            if (r1 == 1){ status[i+j] = true; continue; }
                            if (r1 ==-1){ status[i+j] = 1==line_full(  bls[j], places[j],  place1s[j],  values[j],  0,  false, 0);   continue;  }
                            if (r1 ==-2){   status[i+j] =   1==loop_failure(  bls[j], places[j],  place1s[j],  values[j],  0,  false,  v, 0); continue;   }
                        }
                    }
                    
                    entries_to_lefts[j]  = (int)__builtin_popcountll( (fcas[j][steps[j] ])&m0s[j] ) + steps[j]  * (int)__builtin_popcountll(fcas[j][0]) ;
                    tot_entries[j] = (int)__builtin_popcountll( fcas[j][0]) + (int)__builtin_popcountll(fcas[j][1]  & mask_for_slot_array) ;
                    
                    if (__builtin_expect( (tot_entries[j] > tot_entries_thres ),0)){  // optimize presuming space
                        place1s[j] = h2(places[j],values[j]);
                        status[i+j] = 1==line_full(  bls[j] , places[j] ,  place1s[j] ,  values[j],  0,  false, 0); continue;
                    } else {
                        
                        if (bs[j][start_of_fingerprint_storage+entries_to_lefts[j]] == values[j])  { status[i+j] =  true; continue; }
                        if ( tot_entries[j] -entries_to_lefts[j] >0)
                            memmove( bs[j] + start_of_fingerprint_storage + 1 +  entries_to_lefts[j] ,  bs[j] + start_of_fingerprint_storage + entries_to_lefts[j] , tot_entries[j] -entries_to_lefts[j] );
                        PCF_set1(bs[j],  entries_to_lefts[j]   , values[j]);
                        fcas[j] [steps[j] ]  |=    bit_to_sets[j]  ;
                        status[i+j] =  true;
                    }
                }
            }
        }
    }
    
    inline    int hput_l(int bl,int place, int place1, int filter, uint64_t key, bool early, int depth){
        if( depth > 450){
            return -1;
        }
        int r = CF2_set(bl,place, filter);
        if (r ==1 ) return 1;
        if ( r == 0 ){
            r = CF2_set(bl,place1, filter);
            if (r ==1 ) return 1;
        }
        if ( r == -1){  // we had no space.  Otherwise PIB would have returned 0.
            return line_full(  bl, place,  place1,  filter,  key,  false, depth);
        }
        
        vector < pair < int ,int > > v;
        int r1 = cf_move(bl,place, filter, early,v);
        if (r1 == 1) return 1;
        if (r1 ==-1){
            return line_full(  bl, place,  place1,  filter,  key,  false, depth);
        }
        if (r1 ==-2){
            return   loop_failure(  bl, place,  place1,  filter,  key,  false, v, depth);
        }
        return -1;
    }
    };





#endif /* LCF_h */
