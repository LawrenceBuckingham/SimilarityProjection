#pragma once

#if defined(DO_OMP_TIMER) && DO_OMP_TIMER && USE_OMP
#define OMP_TIMER_DECLARE(name) double hb_elapsed_##name = 0;
#define OMP_TIMER_START(name) double hb_start_##name = omp_get_wtime();
#define OMP_TIMER_END(name) hb_elapsed_##name += omp_get_wtime() - hb_start_##name;
#define OMP_TIMER(name) hb_elapsed_##name
#else
#define OMP_TIMER_DECLARE(name) 
#define OMP_TIMER_START(name) 
#define OMP_TIMER_END(name) 
#define OMP_TIMER(name) 
#endif
