/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   db.hpp
 * Author: Lawrence
 *
 * Created on 22 September 2017, 12:40 AM
 */

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <cstdio>

#ifdef TRON
#undef TRACE
#undef PR
#undef DB
#define TRACE {fprintf(stderr,"%s:%d\n",__FILE__, __LINE__); fflush(stderr);}
#define PR(symbol,format) {fprintf(stderr,"%s:" #format "\n", #symbol, symbol); fflush(stderr);}
#define DB(stuff) stuff
#else
#undef TRACE
#undef PR
#undef DB
#define TRACE 
#define PR(symbol,format) 
#define DB(stuff) 
#endif


