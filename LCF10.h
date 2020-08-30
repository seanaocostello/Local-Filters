//
//  LCF10.h
//  LCF-10bit-conting-42
//


#ifndef LCF10_h
#define LCF10_h

#include <iostream>

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <inttypes.h>
#include <chrono>
#include <list>
#include <vector>
#include <array>
#include  <cstdlib>
#include <iomanip>

#include <cstring>

#include <emmintrin.h>

#include <random>
#include <cmath>


long short_hit = 0;
long short_hit2 = 0;

long lines_full = 0;
long loop_fail = 0;
long removed = 0;
long  moves =0;

long  line_full_ota = 0;
long  line_full_zero_ota = 0;
long line_full_rand  = 0;
long tot46 =0;
long tot45 = 0;
long tot44 = 0;
long found_empty= 0;
long found_nonempty= 0;

int lengths_full[128];
int depth_count[256];
long cucks =0 ;

long lcucks =0 ;
int bad_line = -1 ;
int alt_bad_line = -1;
uint64_t  bad = -1;


bool print1 = false;
bool print = false;

int cuckoo_fail = 0;

long first_fast_fail = 0;
long second_fast_fail = 0;
long first_look_fail = 0;
long first_look_success = 0;
long second_look_fail = 0;
long second_look_success = 0;
long ota_miss=0;

static constexpr int OTA_byte_start  = 10;
static constexpr int OTA_bit_start  = 6;
static constexpr int FCA_byte_start  = 12;

static constexpr int FCA_bit_start  = 6;

static constexpr int number_slots  = 86;
static constexpr int  number_FCA_slots = 41;

static constexpr int LCF10_batch_size = 128; //256;

static constexpr int max_slot_value = 86;

class LCF10 {
public:
    
    unsigned char * ar1;
    unsigned char * ar_end;
    unsigned long long int  hs2[1024];
    long  num_lines;
    long num_lines1;
    int lg_num_lines;
    
    LCF10(   long size){
        
        long sz = size; // one cache line
        if (posix_memalign((void **)&ar1,64, (sz+1) * 64 ) != 0) {
            std::cout << "error on memalign" << std::endl;
            exit(1);
        }
        bzero(ar1,(sz+1)*64);
        num_lines = sz;
        num_lines1 = num_lines-1;
        lg_num_lines = (int)log2(num_lines);
  
        // hs2 is for lines
        for (int i =0;i < 1024; i++){
            hs2[i] = ( (unsigned  int)hash1( rand() +1 )&(num_lines1  ));
            while ((hs2[i]&(num_lines1 )) == 0 )
                hs2[i] = ( (unsigned  int)hash1( rand() +1 )&(  num_lines1))   ;
        }
        
    
     //   print = false;
        
        
    }
    ~LCF10(){
        free(ar1);
    }
    
    
    int PCF_get1(unsigned char * ar, int n, int f){
        int b_start = n*10 + FCA_bit_start;
        int l_start = b_start&(~15); // mask away lower 4 bits, so we are bit aligned on an int
        unsigned  l = *(unsigned  *)(ar + FCA_byte_start +(l_start>>3));
        int shift = b_start&(15);
        unsigned  l1 = (l>>shift);
        return f == (l1&((1<<10)-1));
    }
    
    int PCF_set1(unsigned char * ar, int n, int v ){
        int b_start = n*10 +FCA_bit_start;
        int l_start = b_start&(~15);
        unsigned  * l = (unsigned  *)(ar + FCA_byte_start + (l_start>>3));
        int shift = b_start&(15);
        unsigned bshift = (1<<10)-1;
        unsigned v_clear =  ~(bshift << shift);
        unsigned v_mask = ( (v&bshift) << shift);
        l[0] = ( (l[0]&v_clear) |v_mask);
         if ((ar -ar1)/64  == bad_line)    std::cout << "PCF_set1 " << n << " ret " <<  v << std::endl;
        return 0;
        }
    
         int PCF_get2(unsigned char * ar, int n){
        int b_start = n*10 + FCA_bit_start;
        int l_start = b_start&(~15); // mask away lower 4 bits, so we are bit aligned on a short
        unsigned  l = *(unsigned  *)(ar + FCA_byte_start +(l_start>>3));
        int shift = b_start&(15);
        unsigned  l1 = (l>>shift);
        if ((ar -ar1)/64  == bad_line && (l1&((1<<10)-1)) != 920)    std::cout << "PCF_get2 " << n << " ret " << (l1&((1<<10)-1)) << " bl " <<  (ar -ar1)/64   << std::endl;
        return (l1&((1<<10)-1));
    }
    

    
    
    inline unsigned int CF1_get(int bl,int buc){
        unsigned char * b = ar1 + (bl * 64);
        unsigned long long * fca = (unsigned long long  *)b ;
        int step = buc>>6;
        int off = buc&63;
        unsigned long long bit_to_set = (1UL<<off);
        unsigned long long m0 = ( bit_to_set -1);
        if (! ((fca[step]) &  bit_to_set)){
            return -1;
        }
        int entries_to_left = (int)__builtin_popcountll( (fca[step])&m0) + step * (int)__builtin_popcountll(fca[0]) ;
        return  PCF_get2(b,entries_to_left);
    }
    
    
    
    inline unsigned int CF1_set(int bl,int buc, int v){
        unsigned char * b = ar1 + (bl * 64);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        int step = buc>>6;
        int off = buc&63;
        unsigned long long bit_to_set = (1UL<<off);
        unsigned long long m0 = ( bit_to_set -1);
        int entries_to_left = (int)__builtin_popcountll( (fca[step])&m0) + step * (int)__builtin_popcountll(fca[0]) ;
        int tot_entries = (int)__builtin_popcountll(fca[1] & 0x3FFFFF)  + (int)__builtin_popcountll(fca[0]) ;
        if (tot_entries > number_FCA_slots -1  ){return -1;}
        if ( !((fca[step]) &  bit_to_set) ){
            for (int i = tot_entries -1; i >= entries_to_left ; i--){
                PCF_set1(b,  i + 1  , PCF_get2(b,i));
            }
        }
        PCF_set1(b,  entries_to_left  , v);
        fca[step]  |=    bit_to_set ;
        return 0;
    }
    
    inline unsigned int CF2_set(int bl,int buc, int v){ // return -1 for full, 1 for success, 0 for already occupied
        unsigned char * b = ar1 + (bl * 64);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        int step = buc>>6;
        int off = buc&63;
        unsigned long long  bit_to_set = (1UL<<off);
        unsigned long long m0 = ( bit_to_set -1);
        if ( ((fca[step]) &  bit_to_set) ){return 0;}
        int entries_to_left = (int)__builtin_popcountll( (fca[step])&m0) + step * (int)__builtin_popcountll(fca[0]) ;
        int tot_entries = (int)__builtin_popcountll(fca[1] & 0x3FFFFF)  + (int)__builtin_popcountll(fca[0]) ;
        if (tot_entries > number_FCA_slots -1 ){
            return -1; }
        for (int i =  tot_entries - 1; i >= entries_to_left ; i--){
           PCF_set1(b,  i + 1  , PCF_get2(b,i));
        }
        PCF_set1(b,  entries_to_left  , v);
        fca[step]  |=    bit_to_set ;
        if (bl == bad_line) std::cout << " CF2_set at " << entries_to_left << " to " << v << " pos " << buc << " reget " <<  PCF_get2(b,entries_to_left) <<  std::endl;
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
    
    
     // This computes the  new block.  It needs to
    
    inline  int lh2(  int pl, unsigned int filter){
        unsigned long long  int fh =  ( hs2[filter]) ;
        return (int)(fh^pl) ;
    }
    
    // Warning, if 0x1fff is greater than num_lines, bad things can happen. 8192
    // the in cache involution.  0..95
     inline  int h2(  int pl, unsigned int filter){
         long one = 1;
         long long  int   offset = ((hs2[filter] & 0x3f) +
                                    1) | one; // The offset is at most 63.
         // offset should always be odd.
         offset = (pl & one) ? offset : -offset;
         int64_t output = pl + offset;
         if(output < 0){
             output += max_slot_value;
          }
         if(output >= max_slot_value){
             output -= max_slot_value;
             
         }
         return (int)output;
     }
    
    // This is to replace lh2 when you want a filter with a non power of two number of lines.
    // Warning, if 0x1fff is greater than num_lines, bad things can happen. 8192
    // makes a new line from an old, arbitrary num_lines.
    inline  int lh4(  int pl, unsigned int filter){
        long one = 1;
        long long  int   offset = ((hash1(filter) & 0x1fff) +
                                   1 ) | one;
        // offset should always be odd.
        offset = (pl & one) ? offset : -offset;
        int64_t output = static_cast<int64_t>(pl) + offset;
        if(output < 0){
            output += num_lines;
         }
        if(static_cast<uint64_t>(output) >= num_lines){
            output -= num_lines;
            
        }
        return (int)output;
    }
    
   
    // makes a new place 0..95 from an old.
    inline  int lh3(  int pl, unsigned int filter){
        return pl;
    }
    
    // When you look up an alternate line.  Never needs to look at OTA bits.
    inline  int agg_hget_l_alt( int bl, int place, int place1, int filter){
        unsigned char * b = ar1 + (bl * 64);
        unsigned long bit_to_set;
        int step = place>>6;
        int off = place&63;
        bit_to_set = (1UL<<off);
        int bot_count =(int)__builtin_popcountll(((unsigned long long *)b)[0]);
        if (!((((unsigned long long *)b)[step]) &  bit_to_set) ){
            //   return 0;
        } else {
            unsigned long m0 = ( bit_to_set -1);
            int entries_to_left = (int)__builtin_popcountll( (((unsigned long long *)b)[step])&m0) + step *  bot_count;
            if ( __builtin_expect( ( PCF_get1(b,entries_to_left, filter)) ,0))  return 1;
        }
        int step1 = place1>>6;
        int off1 = place1&63;
        unsigned long bit_to_set1 = (1UL<<off1);
        if (!((((unsigned long long *)b)[step1]) &  bit_to_set1)){
            return 0; }
        unsigned long m01 = ( bit_to_set1 -1);
        int  entries_to_left1 = (int)__builtin_popcountll( (((unsigned long long *)b)[step1])&m01) + step1 * bot_count ;
        if ( __builtin_expect( (PCF_get1(b,entries_to_left1, filter)),0)) return 1;
        return 0;
        
    }
      
    
    
    // The basic lookup, aggressive as it does not use CF1 etc.
    inline  int agg_hget(uint64_t key){
        unsigned long  int h = (unsigned long int) hash1(key);
        int filter = ( h & 1023);
        int place =  (h>>10)&255;
        place = (place*number_slots)>>8;
        int bl = (int)((h>>18)&(num_lines1));
        unsigned char * b = ar1 + (bl * 64);
        int step = place>>6;
        int  off = place&63;
        unsigned long  bit_to_set = (1UL<<off);
       // unsigned short OTA = *(unsigned short *)(b + 10);
        unsigned short OTA = ((*(unsigned *)(b + OTA_byte_start))>>OTA_bit_start)&((1<<16) -1);
        if (OTA == 0){ // We can return 0 out of this block, as the OTA is 0.
            if (!((((unsigned long long *)b)[step]) &  bit_to_set) ){return 0; }
            unsigned long m0 = ( bit_to_set -1);
            int bot_count =(int)__builtin_popcountll(((unsigned long long *)b)[0]);
            int entries_to_left = (int)__builtin_popcountll( (((unsigned long long *)b)[step])&m0) + step *  bot_count;
            if (PCF_get1(b,entries_to_left,filter)) return 1;
            int place1 = h2(place,filter);
            int step1 = place1>>6;
            int off1 = place1&63;
            unsigned long bit_to_set1 = (1UL<<off1);
            if (! ((((unsigned long long *)b)[step1]) &  bit_to_set1)){  return 0; }
            unsigned long m1 = ( bit_to_set1 -1);
            int entries_to_left1 = (int)__builtin_popcountll( (((unsigned long long *)b)[step1])&m1) + step1 * bot_count ;
            if (PCF_get1(b,entries_to_left1, filter)) return 1;
            return 0;
        }
        int place1 = h2(place,filter);
        if (  ((((unsigned long long *)b)[step]) &  bit_to_set) ){
            unsigned long m0;
            m0 = ( bit_to_set -1);
            int entries_to_left = (int)__builtin_popcountll( (((unsigned long long *)b)[step])&m0) + step * (int)__builtin_popcountll(((unsigned long long *)b)[0]) ;
            if (PCF_get1(b,entries_to_left, filter)) return 1;
        }
        off = place1;
        step = off>>6;
        off = off&63;
        bit_to_set = (1UL<<off);
        if ( ((((unsigned long long *)b)[step]) &  bit_to_set)  ){
            unsigned long m0 = ( bit_to_set -1);
            int entries_to_left = (int)__builtin_popcountll( (((unsigned long long *)b)[step])&m0) + step * (int)__builtin_popcountll(((unsigned long long *)b)[0]) ;
            if (  PCF_get1(b,entries_to_left, filter)) return 1;
        }
         if  (1&(OTA>>(filter%16) )){
            int lower = (place<place1)?place:place1;
            int new_bl = lh2(bl,filter);
            int new_place = lh3(lower,filter);
            int new_place1 = h2(new_place, filter);
            return agg_hget_l_alt( new_bl,  new_place, new_place1,  filter);
        }
        return 0;
    }
    
    
    inline  int hget(uint64_t key){
        unsigned long  int h = (unsigned long int) hash1(key);
        int filter = ( h & 1023);
        h =h>>10;
        int place = h&255;
         place = (place*number_slots)>>8;
         h=(h>>8);
        int bl = (int)(h&(num_lines -1 ));
        int place1 = h2(place,filter);
        return hget_l(bl, place, place1, filter, false);
    }
    
    inline  int hget_l( int bl, int place, int place1, int filter , bool secondary){
        int r = CF1_get(bl,place);
        if ( (r&0) == (filter &0)) short_hit++;
        if (r == filter){
            return 1;
        }
        unsigned short OTA = ((*(unsigned *)(ar1 + (bl * 64) + OTA_byte_start))>>OTA_bit_start)&((1<<16) -1);
        if ( r == -1 &&   !secondary
            && OTA == 0 ){
            return 0;
        }
        int r2 = CF1_get(bl,place1 );
        if ( (r2&0) == (filter&0)) short_hit2++;
        if ( r2 == filter){
            return 1;
        }
        //int lower = (place<place1)?place:place1;
        if  ( !secondary && 1&(OTA>>(filter%16) )){ // GGG
            int lower = (place<place1)?place:place1;
            int new_bl = lh2(bl,filter);
            int new_place = lh3(lower,filter);
            int new_place1 = h2(new_place, filter);
            return hget_l( new_bl,  new_place, new_place1,  filter ,true);
        }
        return 0;
    }
    
    
    
    
    
    inline  int  is_room(int bl){
        unsigned char * b = ar1 + (bl * 64);
        unsigned char * fca = b ;
        unsigned long long v0 = *(unsigned long long *)fca;
        unsigned long long v1 = *(unsigned long long *)(fca+8);
        int one_bits_set = (int)__builtin_popcountll(v0);
        int two_bits_set = (int)__builtin_popcountll(v1 & 0x3FFFFF);
        int tot_entries = one_bits_set + two_bits_set ;
        return  42 - tot_entries ;
    }
    
    bool quick_move(int bl, int place, int filt,  std::vector< std::pair < int, int > > &to_move_list ){
        
        int tries = 0;
        to_move_list.push_back(std::make_pair(place, filt));
        int pos = place;
        int filter = CF1_get(bl,pos );
        if (filter == -1){
            return true;
        }
        while(tries++ < 16)
        {
            moves++;
            int next = h2(pos, filter );
            int new_filter =  CF1_get(bl,next ) ;
            if (new_filter == filter) {
                return true;
                break;
            }
            if ( next == place){
                //    int ofil = -1;
                for (auto k: to_move_list){
                    if (place == k.first && filter == k.second){
                        return true;
                    }
                }
                to_move_list.push_back(std::make_pair(next, filter));
                return false;
                
            }
            to_move_list.push_back(std::make_pair(next, filter));
            if(new_filter == -1){  // This needs one more space, so may fail.
                if (is_room(bl) > 0 ) { return true;}
                return false;
            }
            filter = new_filter;
            pos = next;
        }
        return false;  // the fall out case.
    }
    
    inline   int cf_move(int bl, int place, int filt, bool early, std::vector < std::pair < int ,int > > &v_out){
        cucks++;
        if (bl == bad_line)  std::cout << " cf move " << std::endl;
        std::vector< std::pair < int, int > > to_move_list;
        bool rdone = quick_move( bl,  place,  filt, to_move_list );
        if (!rdone) {
            to_move_list.clear();
            rdone = quick_move( bl,  h2(place,filt),  filt, to_move_list );
        }
        if ( rdone){
            for (auto y : to_move_list){
                int r =   CF1_set(bl,y.first,y.second);
                if (r == -1){
                    return -1;
                }
            }
            return 1;
        }
        v_out =  to_move_list;
        cuckoo_fail++;
        return -2;
    }
    
    
    inline  int hdelete(uint64_t key){
        unsigned long  int h  = (unsigned long int) hash1(key);
        int filter = ( h & 1023);
        h = (h>> 10);
        int place = (h&255);  // How do we do this?  Cheap
        place = (place*number_slots)>>8;
        h = (h>>8);
        int bl = (int)(h&(num_lines1));
        int place1 = h2(place,filter);
        return hdelete_l( bl, place,  place1,  filter,  false);
    }
       
       
       inline    int hdelete_l(int bl,int place, int place1, int filter,  bool sec){
           int off = place;
           unsigned char * b = ar1 + (bl * 64);
           unsigned long long *  fca = ( unsigned long long  *)b ;
           unsigned long long v0 = *fca;
           unsigned long long v1 = *(fca+1);
           unsigned long m0;
           unsigned long bit_to_set;
           int step = off>>6;
           off = off&63;
           bit_to_set = (1UL<<off);
           m0 = ( bit_to_set -1);
           if (__builtin_expect( ((fca[step]) &  bit_to_set),0 )){  // already occupied, try alt position optimize for space available.
               int bit_count = (int)__builtin_popcountll(v0);
               int entries_to_left = (int)__builtin_popcountll( (fca[step])&m0) + step * bit_count ;
               int tot_entries = (int)__builtin_popcountll(v0) + (int)__builtin_popcountll(v1 & 0xffffffff) ;
               if (b[14+entries_to_left] == filter) {
                   for (int i = entries_to_left ; i <tot_entries -1; i++){
                       PCF_set1(b,  i  , PCF_get2(b,i+1));
                   }
                   fca[step]  &=    ~bit_to_set ;
                   return 1;
               }
           }
           off = place1;
           step = off>>6;
           off = off&63;
           bit_to_set = (1UL<<off);
           m0 = ( bit_to_set -1);
           if ( __builtin_expect(  ((fca[step]) &  bit_to_set),0 )){ // optimize for space available.
               int bit_count = (int)__builtin_popcountll(v0);
               int entries_to_left = (int)__builtin_popcountll( (fca[step])&m0) + step * bit_count ;
               int tot_entries = (int)__builtin_popcountll(v0) + (int)__builtin_popcountll(v1 & 0xffffffff) ;
               if (b[14+entries_to_left] == filter){
                   for (int i = entries_to_left ; i <tot_entries -1; i++){
                       PCF_set1(b,  i  , PCF_get2(b,i+1));
                   }
                   fca[step]  &=  ~bit_to_set ;
                   return 1;
               }
           }
           if(sec) return 0;
           int lower = (place<place1)?place:place1;
           int new_bl = lh2(bl,filter);
           int new_place = lh3(lower,filter);
           int new_place1 = h2(new_place, filter);
           return hdelete_l( new_bl, new_place,  new_place1,  filter,  true);
           
       }
    
    inline  int hput(uint64_t key){
        unsigned long  int h  = (unsigned long int) hash1(key);
        int filter = ( h & 1023);
        h = (h>> 10);
        int place = (h&255);  // How do we do this?  Cheap
        place = (place*number_slots)>>8;
        h = (h>>8);
        int bl = (int)(h&(num_lines1));
        int place1 = h2(place,filter);
        return  hput_l(bl,place,place1, filter, key, false,0);
        
    }
    
  
    
    int delete_entry( int bl,int place, int & filter){
        removed++;
        int off = place;
        unsigned char * b = ar1 + (bl * 64);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        unsigned long long v0 = *fca;
        unsigned long long v1 = *(fca+1);
        unsigned long long m0;
        unsigned long long bit_to_set;
        int step = off>>6;
        off = off&63;
        bit_to_set = (1UL<<off);
        m0 = ( bit_to_set -1);
        int entries_to_left = (int)__builtin_popcountll( (fca[step])&m0) + step * (int)__builtin_popcountll(v0) ;
        int one_bits_set = (int)__builtin_popcountll(v0);
        int two_bits_set = (int)__builtin_popcountll(v1 & 0x3FFFFF);
        int tot_entries = one_bits_set + two_bits_set ;
        filter = PCF_get2(b,  entries_to_left );
        for (int i = entries_to_left ; i <tot_entries -1; i++){
                if ( i+1 >= 41) std::cout << "delete i = " << i << std::endl;
                        PCF_set1(b,  i  , PCF_get2(b,i+1));
                   }
        fca[step]  &=    ~bit_to_set ;
        lcucks++;
     return 0;
    }
    
    
    inline int is_present(unsigned char * b, int off){
        unsigned long long *  fca = ( unsigned long long  *)b ;
        unsigned long bit_to_set;
        int step = off/64 ; //>>6;
        off = off&63;
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
    
    
    
     int line_full( int bl,int place, int place1, int filter, uint64_t key, bool early, int depth){
         if( depth > 450){
               std::cout << "Depth too high " << depth << " bl " << bl << std::endl;
             return -1;
         }
         lines_full++;
        unsigned char * b = ar1 + (bl * 64);
        unsigned short OTA = ((*(unsigned *)(b + OTA_byte_start))>>OTA_bit_start)&((1<<16) -1);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        int spot = -1;
         int spot1 = -1;
         int first_empty = -1;
         bool emp = false;
         if ( OTA!= 0 && ( depth < 100 || depth%16 == 0)
             ){
              line_full_ota++;
             for (int k = 0; k < number_FCA_slots; k++){
                 int f = PCF_get2(b,k);
                 if (  ((short)(OTA>>(f%16))&1) != 0 ){
                     int new_bl = lh2(bl, PCF_get2(b,k)); //
                     if (PCF_get2(ar1 + new_bl * 64, number_FCA_slots-1 ) == 0){
                        emp = true;
                         found_empty++;
                         int n0 = (int)__builtin_popcountll(fca[0]);
                         spot = ( k >= n0)?(64 + (int)bitselect(fca[1]&0x3FFFFF, k-n0)):(int)bitselect(fca[0], k);
                         break;
                     }
                     
                 } else {
                     if ( first_empty == -1 && k > depth/8 ){
                         int new_bl = lh2(bl, PCF_get2(b,k));//
                         
                       if (PCF_get2(ar1 + new_bl * 64, number_FCA_slots-1 ) == 0){  // good enough to test for empty
                             first_empty = k;
                             int n0 = (int)__builtin_popcountll(fca[0]);
                             spot1 = ( k >= n0)?(64 + (int)bitselect(fca[1]&0x3FFFFF, k-n0)):(int)bitselect(fca[0], k);
                         }
                     }
                 }
             }
         }
         if (!emp && spot != -1) found_nonempty++;
         if (spot == -1 || spot > number_slots-1){
             if ( OTA!= 0 && depth < 16  ) line_full_zero_ota++;
         //    cout << "Chose spot1 " << spot1 << endl;
             spot = spot1;
            }
         if (spot == -1  || spot > number_slots-1){
             //   cout << "Going random as  spot is " << spot  << endl;
               line_full_rand++;
             int buc = rand()%number_slots;
             while(buc > 0 &&  !is_present(b,buc) ){ buc--;  }
             if (buc == 0){
                 buc = number_slots-1;
                 while(buc > 0 &&  !is_present(b,buc) ){ buc --; }
             }
             spot = buc;
             
         }
         int pl = spot;
         int fl1;
         delete_entry(bl,pl, fl1);
         unsigned  * sOTA = (unsigned *)(b + OTA_byte_start);
         hput_l( bl, place,  place1,  filter,  key,  false, depth+1);
         int pl1 = h2(pl,fl1);
         int r =   hput_l_alt( bl, pl,  pl1,  fl1,  0,  false, depth+1);
         if (r == -1) {
             *sOTA |= (1<<(OTA_bit_start + (fl1%16)));
             return -1;
         }
         *sOTA |= (1<<(OTA_bit_start + (fl1%16)));
         return 1;
         
     }
    
     int loop_failure( int bl,int place, int place1, int filter, uint64_t key, bool early, std::vector < std::pair< int, int> > v, int depth){
        loop_fail++;
        int item_no = rand()%(v.size());
        int pl = v[item_no].first;
        int fl1;
         delete_entry(bl,pl, fl1);
        int r =  hput_l( bl, place,  place1,  filter,  key,  false, depth+1);
        unsigned  * sOTA = ( unsigned *)(ar1 + (bl * 64) + OTA_byte_start);
        if (r == -1) {
            *sOTA |= (1<<(6 + (fl1%16)));
            return -1;
        }
        int pl1 =  h2(pl,fl1);
        r = hput_l_alt( bl, pl, pl1,  fl1,  0,  false, depth+1);
        if (r == -1){
              *sOTA |= (1<<(OTA_bit_start + (fl1%16)));
            return -1;
        }
          *sOTA |= (1<<(OTA_bit_start + (fl1%16)));
        return 1;
        
     }
    inline  int hput_l_alt(int bl, int place, int place1, int filter, int key, bool early, int depth){
        int lower = (place<place1)?place:place1;
        int new_bl = lh2(bl,filter);
        int new_place = lh3(lower,filter);
        int new_place1 = h2(new_place, filter);
        if( depth < 256){ depth_count[depth]++; }
        return hput_l( new_bl,  new_place, new_place1,  filter ,key,true, depth+1);  // why true?  Because this is is a secondary insert?  Maybe this is not used.
        
    }

    inline void likely_contains_many( std::vector<uint64_t >& keys,
                                     std::vector<bool>& status,  uint64_t num_keys) {
        for(int i = 0; i < num_keys; i += LCF10_batch_size){
            std::array<uint64_t, LCF10_batch_size> hs;
            std::array<uint32_t, LCF10_batch_size> bls;
            std::array<uint32_t, LCF10_batch_size> values;
            std::array<int, LCF10_batch_size> places;
            std::array<int, LCF10_batch_size> place1s;
            std::array<unsigned char *, LCF10_batch_size> bs;
          //  std::array<uint64_t *, batch_size> fcas;
            std::array<int, LCF10_batch_size> offs;
            std::array<uint64_t, LCF10_batch_size> m0s;
            std::array<uint64_t, LCF10_batch_size> bit_to_sets;
            std::array<int, LCF10_batch_size> steps;
            std::array<int, LCF10_batch_size> entries_to_lefts;
            std::array<bool, LCF10_batch_size> over;
            std::array<short, LCF10_batch_size> OTAs;
            for(int j = 0; j < LCF10_batch_size; j++){
                hs[j] = hash1(keys[i + j]);
            }
            for(int j = 0; j < LCF10_batch_size; j++){
                values[j] = (hs[j])&1023;
                bls[j] = (uint32_t)((hs[j]>>18)&(num_lines1));
                places[j] = ((hs[j]>>10)&255);
                places[j] = (places[j]*max_slot_value)>>8;
            }
            for(int j = 0; j < LCF10_batch_size; j++){
                bs[j] = ar1 + (bls[j] * 64);
                _mm_prefetch((char*)( bs[j] ),_MM_HINT_T0);
            }
            for(int j = 0; j < LCF10_batch_size; j++){
                steps[j] = places[j]>>6;
                offs[j] = places[j]&63;
                bit_to_sets[j] = (1UL<<offs[j]);
            }
            for(int j = 0; j < LCF10_batch_size; j++){
                OTAs[j] = ((*(unsigned *)(bs[j] + OTA_byte_start))>>OTA_bit_start)&((1<<16) -1);
                over[j] = false;
            }
            for(int j = 0; j < LCF10_batch_size; j++){
                if ( OTAs[j] == 0){
                    if (!((((unsigned long long *)bs[j])[steps[j]]) &  bit_to_sets[j]) ){
                        first_fast_fail++; status[i+j] = false;   over[j] = true; }
                    else{
                        m0s[j] = ( bit_to_sets[j] -1);
                        entries_to_lefts[j] = (int)__builtin_popcountll( (((unsigned long long *)bs[j])[steps[j]])&m0s[j]) +
                        steps[j] *  (int)__builtin_popcountll(((unsigned long long *)bs[j])[0]);
                        if ( PCF_get1(bs[j],entries_to_lefts[j],values[j])){
                            first_look_success++; status[i+j] = true;   over[j] = true; }
                        else {
                            first_look_fail++;
                            place1s[j] = h2(places[j],values[j]);
                            offs[j] = place1s[j];
                            steps[j] = offs[j]>>6;
                            offs[j] = offs[j]&63;
                            bit_to_sets[j] = (1UL<<offs[j]);
                            if (! ((((unsigned long long *)bs[j])[steps[j]]) &  bit_to_sets[j])){
                                second_fast_fail++; status[i+j] = false;   over[j] = true; }
                            else {
                                m0s[j] = ( bit_to_sets[j] -1);
                                entries_to_lefts[j] = (int)__builtin_popcountll( (((unsigned long long *)bs[j])[steps[j]])&m0s[j]) +
                                steps[j] *   (int)__builtin_popcountll(((unsigned long long *)bs[j])[0]);
                                status[i+j] = ( PCF_get1(bs[j],entries_to_lefts[j],values[j]));
                                if ( PCF_get1(bs[j],entries_to_lefts[j],values[j])){
                                    second_look_success++; }
                                else { second_look_fail++;}
                                over[j] = true;
                            }
                        }
                    }
                }
            }
            
            for(int j = 0; j < LCF10_batch_size; j++){
                if (!over[j]){
                    place1s[j]  = h2(places[j] ,values[j]);
                    if (  ((((unsigned long long *)bs[j] )[steps[j] ]) &  bit_to_sets[j] ) ){
                        unsigned long m0s[j] ;
                        m0s[j]  = ( bit_to_sets[j]  -1);
                        entries_to_lefts[j]  = (int)__builtin_popcountll( (((unsigned long long *)bs[j] )[steps[j] ])&m0s[j] ) + steps[j]  * (int)__builtin_popcountll(((unsigned long long *)bs[j] )[0]) ;
                        if (PCF_get1(bs[j] ,entries_to_lefts[j],values[j] )) { status[i+j] = true;   over[j] = true; }
                    }
                    if (!over[j]){
                        offs[j]  = place1s[j] ;
                        steps[j]  = offs[j] >>6;
                        offs[j]  = offs[j] &63;
                        bit_to_sets[j]  = (1UL<<offs[j] );
                        if ( ((((unsigned long long *)bs[j] )[steps[j] ]) &  bit_to_sets[j] )  ){
                            m0s[j]  = ( bit_to_sets[j]  -1);
                            entries_to_lefts[j]  = (int)__builtin_popcountll( (((unsigned long long *)bs[j] )[steps[j] ])&m0s[j] ) + steps[j]  * (int)__builtin_popcountll(((unsigned long long *)bs[j] )[0]) ;
                            if (PCF_get1(bs[j] ,entries_to_lefts[j], values[j] )) {  status[i+j] = true;   over[j] = true; }
                        }
                    }
                    if  (!over[j] && (1&(OTAs[j]>>(values[j]%16) ))){ // GGG
                         int lower = (places[j] <place1s[j] )?places[j]:place1s[j] ;
                        int new_bl = lh2(bls[j],values[j]);
                        int new_place = lh3(lower,values[j]);
                        int new_place1 = h2(new_place, values[j]);
                         status[i+j] = hget_l( new_bl,  new_place, new_place1,  values[j] ,true);
                    } else if  (!over[j] ){
                        ota_miss++;
                        status[i+j] = false;
                    }
                }
            }
        }
    }
    
    
    inline  int agg_hput(uint64_t key){
        unsigned long  int h  = (unsigned long int) hash1(key);
        int filter = ( h & 1023);
        h = (h>> 10);
        int place = (h&255);
        place = (place*number_slots)>>8;
        h = (h>>8);
        int bl = (int)(h&(num_lines1));
        int place1 = h2(place,filter);
        return  agg_hput_l(bl,place,place1, filter, key, false,0);
        
    }
    
    
    inline    int agg_hput_l(int bl,int place, int place1, int filter, uint64_t key, bool early, int depth){
        int off = place;
        unsigned char * b = ar1 + (bl * 64);
        unsigned long long *  fca = ( unsigned long long  *)b ;
        unsigned long long v0 = *fca;
        unsigned long long v1 = *(fca+1);
        unsigned long m0;
        unsigned long bit_to_set;
        int step = off>>6;
        off = off&63;
        bit_to_set = (1UL<<off);
        m0 = ( bit_to_set -1);
        if (__builtin_expect( ((fca[step]) &  bit_to_set),0 )){  // already occupied, try alt position optimize for space available.
            int bit_count = (int)__builtin_popcountll(v0);
            int entries_to_left = (int)__builtin_popcountll( (fca[step])&m0) + step * bit_count ;
            if ( PCF_get1(b,entries_to_left, filter)) return 1;
            off = place1;
            step = off>>6;
            off = off&63;
            bit_to_set = (1UL<<off);
            m0 = ( bit_to_set -1);
            if ( __builtin_expect(  ((fca[step]) &  bit_to_set),0 )){ // optimize for space available.
                int entries_to_left = (int)__builtin_popcountll( (fca[step])&m0) + step * bit_count ;
                 if ( PCF_get1(b,entries_to_left, filter)) return 1;
                int tot_entries = (int)__builtin_popcountll(v0) + (int)__builtin_popcountll(v1 & 0x3FFFFF) ;
                if (__builtin_expect(  (tot_entries > number_FCA_slots -1 ),0)){ return line_full(  bl, place,  place1,  filter,  key,  false, depth); } // optimize presuming space
                std::vector < std::pair < int ,int > > v;
                int r1 = cf_move(bl,place, filter, early,v);
                if (r1 == 1) return 1;
                if (r1 ==-1){ return line_full(  bl, place,  place1,  filter,  key,  false, depth);    }
                if (r1 ==-2){   return   loop_failure(  bl, place,  place1,  filter,  key,  false,  v, depth);   }
                return -1;
            }
        }
        int entries_to_left = (int)__builtin_popcountll( (fca[step])&m0) + step * (int)__builtin_popcountll(v0) ;
        int one_bits_set = (int)__builtin_popcountll(v0);
        int two_bits_set = (int)__builtin_popcountll(v1 & 0x3FFFFF);
        int tot_entries = one_bits_set + two_bits_set ;
        if (__builtin_expect( (tot_entries > number_FCA_slots -1  ),0)){  // optimize presuming space
            return line_full(  bl, place,  place1,  filter,  key,  false, depth);
        } else {
            PCF_set1(b,  entries_to_left  , filter);
            fca[step]  |=    bit_to_set ;
            return 1;
        }
    }
    
    
    inline void add_many( std::vector<uint64_t >& keys,
                         std::vector<bool>& status,  uint64_t num_keys) {
        for(int i = 0; i < num_keys; i += LCF10_batch_size){
            std::array<uint64_t, LCF10_batch_size> hs;
            std::array<int, LCF10_batch_size> bls;
            std::array<uint32_t, LCF10_batch_size> values;
            std::array<int, LCF10_batch_size> places;
            std::array<int, LCF10_batch_size> place1s;
            std::array<unsigned char *, LCF10_batch_size> bs;
            std::array<uint64_t *, LCF10_batch_size> fcas;
            std::array<int, LCF10_batch_size> offs;
            std::array<uint64_t, LCF10_batch_size> m0s;
            std::array<uint64_t, LCF10_batch_size> bit_to_sets;
            std::array<int, LCF10_batch_size> steps;
            std::array<int, LCF10_batch_size> entries_to_lefts;
            std::array<int, LCF10_batch_size> tot_entries;
            std::array<bool, LCF10_batch_size> over;
            std::array<uint64_t, LCF10_batch_size> pres;
            std::array<bool, LCF10_batch_size> easy;
            for(int j = 0; j < LCF10_batch_size; j++){
                hs[j] = hash1(keys[i + j]);
                pres[j] = 0;
                easy[j] = false;
            }
            for(int j = 0; j < LCF10_batch_size; j++){
                values[j] = (hs[j])&255;
                bls[j] = (int)((hs[j]>>17)&(num_lines1));
                places[j] = ((hs[j]>>10)&255);
                places[j] = (places[j]*max_slot_value)>>8;
                offs[j] = places[j];
                bs[j] = ar1 + (bls[j] * 64);
                fcas[j] = ( uint64_t *)bs[j] ;
                steps[j] = (offs[j]/64);
                offs[j] = offs[j]&63;
                bit_to_sets[j] = (1UL<<offs[j]);
                m0s[j] = ( bit_to_sets[j] -1);
                over[j] = false;
                entries_to_lefts[j]  = (int)__builtin_popcountll( (fcas[j][steps[j] ])&m0s[j] ) + steps[j]  * (int)__builtin_popcountll(fcas[j][0]) ;
                tot_entries[j] = (int)__builtin_popcountll( fcas[j][0]) + (int)__builtin_popcountll(fcas[j][1] & 0x3FFFFF ) ;
            }
            for(int j = 0; j < LCF10_batch_size; j++){
                if  ( __builtin_expect( ( ( !(1UL&( pres[ 127&(bls[j]>>8) ]>>(bls[j]&63))))
                                         && !  ((fcas[j][steps[j]]) &  bit_to_sets[j]) && tot_entries[j] <= number_FCA_slots -1  )   ,1)){
                    easy[j] = true;
                    pres[ 127&(bls[j]>>8)] |= (1UL<<(bls[j]&63));
                    if (__builtin_expect( PCF_get1(bs[j],entries_to_lefts[j],values[j]),0)) { status[i+j] =  true ;continue; }
                    for (int k = tot_entries[j] -1; k >= entries_to_lefts[j] ; k--){
                        PCF_set1(bs[j],  k + 1  , PCF_get2(bs[j],k));
                    }
                    PCF_set1(bs[j],  entries_to_lefts[j] , values[j]);
                    fcas[j] [steps[j] ]  |=    bit_to_sets[j]  ;
                    status[i+j] = true;
                }
            }
            for(int j = 0; j < LCF10_batch_size; j++){
                if (__builtin_expect( (!easy[j]),0)){
                    place1s[j] = h2(places[j],values[j]);
                    status[i+j] = (1 == hput_l(bls[j], places[j],place1s[j], values[j],0,false, 0));
                    continue;
                    if (__builtin_expect( ((fcas[j][steps[j]]) &  bit_to_sets[j]),0 )){  // already occupied, try alt position optimize for space available.
                        tot_entries[j] = (int)__builtin_popcountll( fcas[j][0]) + (int)__builtin_popcountll(fcas[j][1] & 0x3FFFFF ) ;
                        entries_to_lefts[j] = (int)__builtin_popcountll( (fcas[j][steps[j]])&m0s[j]) + steps[j] * (int)__builtin_popcountll(fcas[j][0]) ;
                        if (PCF_get1(bs[j],entries_to_lefts[j],values[j])) { status[i+j] =  true;continue; }
                        place1s[j] = h2(places[j],values[j]);
                        steps[j] = (offs[j]/64);
                        offs[j] = place1s[j]&63;
                        bit_to_sets[j] = (1UL<<offs[j]);
                        m0s[j] = ( bit_to_sets[j] -1);
                        if ( __builtin_expect(  ((fcas[j][steps[j]]) &  bit_to_sets[j]),0 )){ // optimize for space available.
                            entries_to_lefts[j] = (int)__builtin_popcountll( (fcas[j][steps[j]])&m0s[j]) + steps[j] * (int)__builtin_popcountll(fcas[j][0]) ;
                            if (PCF_get1(bs[j],entries_to_lefts[j],values[j]))  { status[i+j] =  true; continue; }
                            tot_entries[j] = (int)__builtin_popcountll( fcas[j][0]) + (int)__builtin_popcountll(fcas[j][1] & 0x3FFFFF) ;
                            if (__builtin_expect(  (tot_entries[j] > number_FCA_slots -1 ),0)){ status[i+j] =  1==line_full(  bls[j], places[j],  place1s[j],  values[j],  0,  false, 0); } // optimize presuming space
                            std::vector < std::pair < int ,int > > v;
                            int r1 = cf_move(bls[j],places[j], values[j], false,v);
                            if (r1 == 1){ status[i+j] = true; continue; }
                            if (r1 ==-1){ status[i+j] = 1==line_full(  bls[j], places[j],  place1s[j],  values[j],  0,  false, 0);   continue;  }
                            if (r1 ==-2){   status[i+j] =   1==loop_failure(  bls[j], places[j],  place1s[j],  values[j],  0,  false,  v, 0); continue;   }
                        }
                    }
                    entries_to_lefts[j]  = (int)__builtin_popcountll( (fcas[j][steps[j] ])&m0s[j] ) + steps[j]  * (int)__builtin_popcountll(fcas[j][0]) ;
                    tot_entries[j] = (int)__builtin_popcountll( fcas[j][0]) + (int)__builtin_popcountll(fcas[j][1]  & 0x3FFFFF) ;
                    if (__builtin_expect( (tot_entries[j] > number_FCA_slots-1  ),0)){  // optimize presuming space
                        place1s[j] = h2(places[j],values[j]);
                        status[i+j] = 1==line_full(  bls[j] , places[j] ,  place1s[j] ,  values[j],  0,  false, 0); continue;
                    } else {
                        if (PCF_get1(bs[j],entries_to_lefts[j],values[j]))  { status[i+j] =  true; continue; }
                        if ( tot_entries[j] -entries_to_lefts[j] >0)
                              for (int k = tot_entries[j] -1; k >= entries_to_lefts[j] ; k--){
                                                  PCF_set1(bs[j],  k + 1  , PCF_get2(bs[j],k));
                                              }
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
              std::cout << "Depth too high " << depth << " bl " << bl << std::endl;
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
        
        std::vector < std::pair < int ,int > > v;
        int r1 = cf_move(bl,place, filter, early,v);
        if (r1 == 1) return 1;
        if (r1 ==-1){
         int r3= line_full(  bl, place,  place1,  filter,  key,  false, depth);
            if (r3 == -1)std::cout << " HPUT_L:  line full  failure full " << std::endl;
            return r3;
        }
        if (r1 ==-2){
            int r2 =   loop_failure(  bl, place,  place1,  filter,  key,  false, v, depth);
            if(r2 ==-1)     std::cout << " HPUT_L:  loop failure full  bl " << bl << std::endl;
            return r2;
            }
        return -1;
    }
    };




#endif /* LCF10_h */
