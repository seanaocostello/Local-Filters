//
//  main.cpp
//  Test Filters
//


#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <inttypes.h>
#include <chrono>
#include <list>
#include <array>
#include  <cstdlib>
#include <iomanip>

#include <algorithm>    // std::sort
#include <vector>       // std::vector

#include <cstring>

#include <emmintrin.h>

#include <random>
#include <cmath>

using namespace std;



#ifndef LL
int repeat_times = 1;


#include "LCF.h"
#include "PQF.h"
#include "LCF10.h"
#else

int repeat_times = 11;

#include "LCF.h"
#include "PQF.h"
#include "LCF10.h"
#endif

bool sec (pair<double, double> i, pair<double, double>  j) { return (i.second > j.second ); }



int main(int argc, const char * argv[]) {
  
    
    
    // Basic idea is
    // For vacuum 12 bit, vacuum-fan 12bit
    // vacuum 13 bit and vacuum fan 13 bit
    // LCF
    // PQF
    // LCF10
    
    // We do the folllowing
    
    std::mt19937_64 generator(0);
    
    int fac = 1024 * 1024 * 4; //
    std::vector<double> lfs = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
        0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.96,
        0.97, 0.98, 0.99};
    
   
        vector<pair<double,double> > run_times_lcf;
         vector<pair<double,double> > run_times_lcf10;
         vector<pair<double,double> > run_times_pqf;
        vector<pair<double,double> > run_times_vf12;
        vector<pair<double,double> > run_times_pvf12;
        vector<pair<double,double> > run_times_vf13;
        vector<pair<double,double> > run_times_pvf13;
        
        vector<pair<double,double> > false_positive_lcf;
        vector<pair<double,double> > false_positive_lcf10;
        vector<pair<double,double> > false_positive_pqf;
        vector<pair<double,double> > false_positive_vf12;
        vector<pair<double,double> > false_positive_pvf12;
        vector<pair<double,double> > false_positive_vf13;
        vector<pair<double,double> > false_positive_pvf13;
    cout << "LCF " << endl;
        
        for(auto i : lfs){
        for (int y = 0; y < 1; y++){
            //fill to a certain percentage
            vector<uint64_t> v,v1,v41;
            long top = 128 * ((int)(46 * fac * i)/128); // 46 for data
            long top_41 = 128 * ((int)(41 * fac * i)/128); // 46 for data
            for (int j =0 ;j < top;j++) v.push_back(generator());
            for (int j =0 ;j < top_41;j++) v41.push_back(generator());
            int test_length = 1024 * 1024; //128 * ((int)(top * 0.001 )/128);
            for (int j =0 ;j <  test_length ;j++) v1.push_back(generator());  // 0.1% of the total, so
        
            
            for (int z = 0; z < repeat_times; z++){
            // First LCF
              LCF * bf = new LCF( fac ); //
            for (auto c : v) bf->hput(c);
            for (auto c : v) bf->hget(c);
            vector<bool> status;
            status.reserve(test_length );
            for (int k = 0; k < test_length; k++){  status[k] = false;}
              bf->likely_contains_many(v1,status, test_length);
            auto   start = std::chrono::high_resolution_clock::now();
            bf->likely_contains_many(v1,status, test_length);
            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            run_times_lcf.push_back(make_pair(i,v1.size()/ (1024 * 1024 * elapsed.count()   )));
          //  cout <<  " ( " <<  i  << " , " << v1.size()/ (1024 * 1024 * elapsed.count() ) <<  " ) " << endl;
          //    cout << " tl " << test_length << " time " << elapsed.count() << endl;
            delete bf;
            }
            std::sort (run_times_lcf.begin(), run_times_lcf.end(), sec);
            cout <<  "  ( " <<  run_times_lcf[0].first   << " , " <<   run_times_lcf[0].second  <<  " ) " << endl;
            
            run_times_lcf.clear();
        }}
       cout << "PQF  " << endl;

            for(auto i : lfs){
            for (int y = 0; y < 1; y++){
                //fill to a certain percentage
                vector<uint64_t> v,v1,v41;
                long top = 128 * ((int)(46 * fac * i)/128); // 46 for data
                long top_41 = 128 * ((int)(41 * fac * i)/128); // 46 for data
                for (int j =0 ;j < top;j++) v.push_back(generator());
                for (int j =0 ;j < top_41;j++) v41.push_back(generator());
                int test_length = 1024 * 1024; //128 * ((int)(top * 0.001 )/128);
                for (int j =0 ;j <  test_length ;j++) v1.push_back(generator());  // 0.1% of the total, so
            
             for (int z = 0; z < repeat_times; z++){
            // PQF
            PQF * bf1 = new PQF( fac ); //
            for (auto c : v) bf1->RS_setc(c);
            vector<bool> status1;
            status1.reserve(test_length );
           auto start = std::chrono::high_resolution_clock::now();
            bf1->likely_contains_many(v1,status1, test_length);
           auto finish = std::chrono::high_resolution_clock::now();
           std::chrono::duration<double>  elapsed = finish - start;
            run_times_pqf.push_back(make_pair(i,v1.size()/ (1024 * 1024 * elapsed.count()   )));
          //  cout <<  " ( " <<  i  << " , " << v1.size()/ (1024 * 1024 * elapsed.count() ) <<  " ) " << endl;
                    delete bf1;
             }
                std::sort (run_times_pqf.begin(), run_times_pqf.end(), sec);
                           cout <<  "   ( " <<  run_times_pqf[0].first   << " , " <<   run_times_pqf[0].second  <<  " ) " << endl;
                run_times_pqf.clear();
                
            }}
      
  
  
       cout << "LCF10 lookup  " << endl;
    for(auto i : lfs){
        for (int y = 0; y < 1; y++){
            //fill to a certain percentage
            vector<uint64_t> v,v1,v41;
            long top = 128 * ((int)(46 * fac * i)/128); // 46 for data
            long top_41 = 128 * ((int)(41 * fac * i)/128); // 46 for data
            for (int j =0 ;j < top;j++) v.push_back(generator());
            for (int j =0 ;j < top_41;j++) v41.push_back(generator());
            int test_length = 1024 * 1024; //128 * ((int)(top * 0.001 )/128);
            for (int j =0 ;j <  test_length ;j++) v1.push_back(generator());  // 0.1% of the total, so
            
            
            for (int z = 0; z < repeat_times; z++){
                //   cout << "Starting LCF10 " << endl;
                // First LCF10
                LCF10 * bf2 = new LCF10( fac ); //
                for (auto c : v41) bf2->hput(c);
                vector<bool> status2;
                status2.reserve(test_length );
                auto start = std::chrono::high_resolution_clock::now();
                bf2->likely_contains_many(v1,status2, test_length);
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double>  elapsed = finish - start;
                //  cout <<  "   ( " <<  i  << " , " << v1.size()/ (1024 * 1024 * elapsed.count() ) <<  " ) " << endl;
                run_times_lcf10.push_back(make_pair(i,v1.size()/ (1024 * 1024 * elapsed.count()   )));
                delete bf2;
            }
            std::sort (run_times_lcf10.begin(), run_times_lcf10.end(), sec);
            cout <<  "   ( " <<  run_times_lcf10[0].first   << " , " <<   run_times_lcf10[0].second  <<  " ) " << endl;
            run_times_lcf10.clear();
        }}
        
    //  I should probably add insert and delete ..
    cout << "LCF insert " << endl;
    run_times_lcf.clear();
    for(auto i : lfs){
        for (int y = 0; y < 1; y++){
            //fill to a certain percentage
            vector<uint64_t> v,v1,v41;
            long top = 128 * ((int)(46 * fac * i)/128); // 46 for data
         //   long top_41 = 128 * ((int)(41 * fac * i)/128); // 46 for data
            for (int j =0 ;j < top;j++) v.push_back(generator());
         //   for (int j =0 ;j < top_41;j++) v41.push_back(generator());
            int test_length = 128 * ((int)(top * 0.001 )/128);
            for (int j =0 ;j <  test_length ;j++) v1.push_back(generator());  // 0.1% of the total, so
            
            
            for (int z = 0; z < repeat_times; z++){
                // First LCF
                LCF * bf = new LCF( fac ); //
                for (auto c : v) bf->hput(c);
                for (auto c : v) bf->hget(c);
                vector<bool> status;
                status.reserve(test_length );
                for (int k = 0; k < test_length; k++){  status[k] = false;}
                // bf->likely_contains_many(v1,status, test_length);
                auto   start = std::chrono::high_resolution_clock::now();
                bf->add_many(v1,status, test_length);
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = finish - start;
                run_times_lcf.push_back(make_pair(i,v1.size()/ (1024 * 1024 * elapsed.count()   )));
                //  cout <<  " ( " <<  i  << " , " << v1.size()/ (1024 * 1024 * elapsed.count() ) <<  " ) " << endl;
                //    cout << " tl " << test_length << " time " << elapsed.count() << endl;
                delete bf;
            }
            std::sort (run_times_lcf.begin(), run_times_lcf.end(), sec);
            cout <<  "  ( " <<  run_times_lcf[0].first   << " , " <<   run_times_lcf[0].second  <<  " ) " << endl;
            
            run_times_lcf.clear();
        }}
    
    
    
    
    cout << "PQF insert " << endl;
    
    for(auto i : lfs){
        for (int y = 0; y < 1; y++){
            //fill to a certain percentage
            vector<uint64_t> v,v1,v41;
            long top = 128 * ((int)(46 * fac * i)/128); // 46 for data
            long top_41 = 128 * ((int)(41 * fac * i)/128); // 46 for data
            for (int j =0 ;j < top;j++) v.push_back(generator());
            for (int j =0 ;j < top_41;j++) v41.push_back(generator());
            int test_length = 128 * ((int)(top * 0.001 )/128);
            for (int j =0 ;j <  test_length ;j++) v1.push_back(generator());  // 0.1% of the total, so
            
            for (int z = 0; z < repeat_times; z++){
                // PQF
                PQF * bf1 = new PQF( fac ); //
                for (auto c : v41) bf1->RS_setc(c);
                vector<int> status1;
                status1.reserve(test_length );
                auto start = std::chrono::high_resolution_clock::now();
                bf1->add_many(v1,status1, test_length);
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double>  elapsed = finish - start;
                run_times_pqf.push_back(make_pair(i,v1.size()/ (1024 * 1024 * elapsed.count()   )));
                //  cout <<  " ( " <<  i  << " , " << v1.size()/ (1024 * 1024 * elapsed.count() ) <<  " ) " << endl;
                delete bf1;
            }
            std::sort (run_times_pqf.begin(), run_times_pqf.end(), sec);
            cout <<  "   ( " <<  run_times_pqf[0].first   << " , " <<   run_times_pqf[0].second  <<  " ) " << endl;
            run_times_pqf.clear();
            
        }}
      
    cout << "LCF10 insert " << endl;
    for(auto i : lfs){
        for (int y = 0; y < 1; y++){
            //fill to a certain percentage
            vector<uint64_t> v,v1,v41;
         //   long top = 128 * ((int)(46 * fac * i)/128); // 46 for data
            long top_41 = 128 * ((int)(41 * fac * i)/128); // 46 for data
           //for (int j =0 ;j < top;j++) v.push_back(generator());
            for (int j =0 ;j < top_41;j++) v41.push_back(generator());
            int test_length = 128 * ((int)(top_41 * 0.001 )/128);
            for (int j =0 ;j <  test_length ;j++) v1.push_back(generator());  // 0.1% of the total, so
            
            
            for (int z = 0; z < repeat_times; z++){
                //   cout << "Starting LCF10 " << endl;
                // First LCF10
                LCF10 * bf2 = new LCF10( fac ); //
                for (auto c : v41) bf2->hput(c);
                vector<bool> status2;
                status2.reserve(test_length );
                auto start = std::chrono::high_resolution_clock::now();
                bf2->add_many(v1,status2, test_length);
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double>  elapsed = finish - start;
                //  cout <<  "   ( " <<  i  << " , " << v1.size()/ (1024 * 1024 * elapsed.count() ) <<  " ) " << endl;
                run_times_lcf10.push_back(make_pair(i,v1.size()/ (1024 * 1024 * elapsed.count()   )));
                delete bf2;
            }
            std::sort (run_times_lcf10.begin(), run_times_lcf10.end(), sec);
            cout <<  "   ( " <<  run_times_lcf10[0].first   << " , " <<   run_times_lcf10[0].second  <<  " ) " << endl;
            run_times_lcf10.clear();
        }}
    
    
    //  I should probably add insert and delete ..
    cout << " LCF delete " << endl;
    run_times_lcf.clear();
    for(auto i : lfs){
        for (int y = 0; y < 1; y++){
            //fill to a certain percentage
            vector<uint64_t> v,v1,v41;
            long top = 128 * ((int)(46 * fac * i)/128); // 46 for data
            long top_41 = 128 * ((int)(41 * fac * i)/128); // 46 for data
            for (int j =0 ;j < top;j++) v.push_back(generator());
            for (int j =0 ;j < top_41;j++) v41.push_back(generator());
            int test_length = 64 * 1024; //128 * ((int)(top * 0.001 )/128);
            for (int j =0 ;j <  test_length ;j++) v1.push_back(generator());  // 0.1% of the total, so
            
            
            for (int z = 0; z < repeat_times; z++){
                // First LCF
                LCF * bf = new LCF( fac ); //
                for (auto c : v) bf->hput(c);
                for (auto c : v) bf->hget(c);
                vector<bool> status;
                status.reserve(test_length );
                for (int k = 0; k < test_length; k++){  status[k] = false;}
                // bf->likely_contains_many(v1,status, test_length);
                auto   start = std::chrono::high_resolution_clock::now();
                for (int k = 0; k < test_length; k++){  bf->hdelete(v1[k]); }
                //bf->add_many(v1,status, test_length);
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = finish - start;
                run_times_lcf.push_back(make_pair(i,v1.size()/ (1024 * 1024 * elapsed.count()   )));
                //  cout <<  " ( " <<  i  << " , " << v1.size()/ (1024 * 1024 * elapsed.count() ) <<  " ) " << endl;
                //    cout << " tl " << test_length << " time " << elapsed.count() << endl;
                delete bf;
            }
            std::sort (run_times_lcf.begin(), run_times_lcf.end(), sec);
            cout <<  "  ( " <<  run_times_lcf[0].first   << " , " <<   run_times_lcf[0].second  <<  " ) " << endl;
            
            run_times_lcf.clear();
        }}
    
    
    
    
    cout << "PQF delete " << endl;
    
    for(auto i : lfs){
        for (int y = 0; y < 1; y++){
            //fill to a certain percentage
            vector<uint64_t> v,v1,v41;
            long top = 128 * ((int)(46 * fac * i)/128); // 46 for data
            long top_41 = 128 * ((int)(41 * fac * i)/128); // 46 for data
            for (int j =0 ;j < top;j++) v.push_back(generator());
            for (int j =0 ;j < top_41;j++) v41.push_back(generator());
            int test_length = 64 * 1024; //128 * ((int)(top * 0.001 )/128);
            for (int j =0 ;j <  test_length ;j++) v1.push_back(generator());  // 0.1% of the total, so
            
            for (int z = 0; z < repeat_times; z++){
                // PQF
                PQF * bf1 = new PQF( fac ); //
                for (auto c : v41) bf1->RS_setc(c);
                vector<int> status1;
                status1.reserve(test_length );
                auto start = std::chrono::high_resolution_clock::now();
                 for (int k = 0; k < test_length; k++){  bf1->RS_delete(v1[k]); }
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double>  elapsed = finish - start;
                run_times_pqf.push_back(make_pair(i,v1.size()/ (1024 * 1024 * elapsed.count()   )));
                //  cout <<  " ( " <<  i  << " , " << v1.size()/ (1024 * 1024 * elapsed.count() ) <<  " ) " << endl;
                delete bf1;
            }
            std::sort (run_times_pqf.begin(), run_times_pqf.end(), sec);
            cout <<  "   ( " <<  run_times_pqf[0].first   << " , " <<   run_times_pqf[0].second  <<  " ) " << endl;
            run_times_pqf.clear();
            
        }}
      
    cout << "LCF10 delete " << endl;
    for(auto i : lfs){
        for (int y = 0; y < 1; y++){
            //fill to a certain percentage
            vector<uint64_t> v,v1,v41;
            long top = 128 * ((int)(46 * fac * i)/128); // 46 for data
            long top_41 = 128 * ((int)(41 * fac * i)/128); // 46 for data
            for (int j =0 ;j < top;j++) v.push_back(generator());
            for (int j =0 ;j < top_41;j++) v41.push_back(generator());
            int test_length = 64 * 1024; //128 * ((int)(top * 0.001 )/128);
            for (int j =0 ;j <  test_length ;j++) v1.push_back(generator());  // 0.1% of the total, so
            
            
            for (int z = 0; z < repeat_times; z++){
                //   cout << "Starting LCF10 " << endl;
                // First LCF10
                LCF10 * bf2 = new LCF10( fac ); //
                for (auto c : v41) bf2->hput(c);
                vector<bool> status2;
                status2.reserve(test_length );
                auto start = std::chrono::high_resolution_clock::now();
                 for (int k = 0; k < test_length; k++){  bf2->hdelete(v1[k]); }
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double>  elapsed = finish - start;
                //  cout <<  "   ( " <<  i  << " , " << v1.size()/ (1024 * 1024 * elapsed.count() ) <<  " ) " << endl;
                run_times_lcf10.push_back(make_pair(i,v1.size()/ (1024 * 1024 * elapsed.count()   )));
                delete bf2;
            }
            std::sort (run_times_lcf10.begin(), run_times_lcf10.end(), sec);
            cout <<  "   ( " <<  run_times_lcf10[0].first   << " , " <<   run_times_lcf10[0].second  <<  " ) " << endl;
            run_times_lcf10.clear();
        }}
    
    
    
    
    
    return 0;
}
