// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2016 AURA/LSST.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#ifndef LSST_MEAS_BASE_libfftw_h_INCLUDED
#define LSST_MEAS_BASE_libfftw_h_INCLUDED

#include "fftw3.h"

#include "lsst/pex/logging.h"

/* 
 * A class returning an object with function pointers to libfftw3 API.
 *
 * Multiple implementations exist of FFTW3 API for Fast Fourier Transforms. 
 * The most well known one is libfftw3 itself, but there are also Intel MKL
 * fftw3 API wrappers.  Unfortunately, the MKL implementation of FFTW3 APIs
 * is incomplete, and returns NULL plans in some cases (e.g., Intel MKL does
 * not support transforms in half-complex format).
 *
 * While we always link our code against libfftw3, it's still possible that
 * at runtime it gets linked against the MKL version.  If that happens, and
 * we attempt using an unimplemented transform, our code will crash (at
 * best) or just return an incorrect result.  As far as we can tell, this is
 * only a problem on Linux with Anaconda, which builds numpy against MKL.
 *
 * This class makes sure that its member function pointers always resolve to
 * the libfftw3 implementations.  If your code requires libfftw3
 * functionality, make FFTW3 calls through an instance of this object.  For
 * example:
 *
 *     auto libfftw3 = libfftw3_t::api();
 *
 *     fftw_plan plan = libfftw3->fftw_plan_r2r_2d(wid, wid, c, c, FFTW_R2HC, FFTW_R2HC, FFTW_ESTIMATE);
 *     libfftw3->fftw_execute(plan);
 *     libfftw3->fftw_destroy_plan(plan);
 *
 * Note that if you create a plan using libfftw3_t functions, you must also
 * use and destroy it with these functions (and vice versa). Otherwise your
 * code is likely to crash.
 *
 */

#define DECLARE_FFTW(func) decltype(&::func) func;

#ifdef __linux__
    #define BIND_FFTW(func)    this->func = (decltype(&::func))dlsym(dlhandle, #func);

    #include <dlfcn.h>
#else
    #define BIND_FFTW(func)    this->func = ::func;
#endif

struct libfftw3_t
{
private:
    void *dlhandle;

public:
    // Functions from libfftw3 that we want to expose.
    //
    // Don't forget to add a corresponding BIND_FFTW() entry to the
    // constructor.
    //
    DECLARE_FFTW(fftw_plan_r2r_1d)
    DECLARE_FFTW(fftw_plan_r2r_2d)
    DECLARE_FFTW(fftw_plan_r2r_3d)
    DECLARE_FFTW(fftw_plan_r2r)
    DECLARE_FFTW(fftw_execute)
    DECLARE_FFTW(fftw_destroy_plan)

public:
    static const libfftw3_t *api()
    {
        // Return the pointer to libfftw3 API. Recommended use:
        //
        //     auto libfftw3 = libfftw3_t::api();
        //

        static std::auto_ptr<libfftw3_t> singleton(NULL);

        if (!singleton.get())
        {
            singleton.reset(new libfftw3_t);
        }

        return singleton.get();
    }

private:
    libfftw3_t() : dlhandle(NULL)
    {
        // Rebind the symbols using the correct library
        dlhandle = get_fftw_lib_handle();

        // Do the actual symbol finding and binding
        BIND_FFTW(fftw_plan_r2r_1d)
        BIND_FFTW(fftw_plan_r2r_2d)
        BIND_FFTW(fftw_plan_r2r_3d)
        BIND_FFTW(fftw_plan_r2r)
        BIND_FFTW(fftw_execute)
        BIND_FFTW(fftw_destroy_plan)
    }

    #ifdef __linux__
    ~libfftw3_t()
    {
        if(dlhandle)
            dlclose(dlhandle);
    }
    #endif

    void *get_fftw_lib_handle()
    {
        /*
            Find (or load) the real libfftw3 and bind its symbols.
        */
        #ifdef __linux__
            namespace pexLog = lsst::pex::logging;
            pexLog::Log log(pexLog::Log::getDefaultLog(), "lsst.meas.base.libfftw3_t");

            // Test if there are MKL symbols in the same library where we found FFTW symbols.
            // If this is the case, then MKL is overriding the real libfftw.
            void* func = dlsym(RTLD_DEFAULT, "fftw_execute");
            Dl_info dlInfo;
            dladdr(func, &dlInfo);
            void *dlhandle = dlopen(dlInfo.dli_fname, RTLD_LOCAL | RTLD_DEEPBIND | RTLD_LAZY);
            if(dlsym(dlhandle, "mkl_dft_dfticomputeforward")) {
                // MKL has overridden FFTW. Let's find our real FFTW library and
                // load it instead.
                dlclose(dlhandle);

                // Find the shadowed function
                void* func = dlsym(RTLD_NEXT, "fftw_execute");
                dladdr(func, &dlInfo);

                // Open libfftw3 with dlmopen, LM_ID_NEWLM, and
                // RTLD_DEEPBIND, ensuring we get all the original libfftw3
                // functions and nothing leaks in from the already-loaded MKL.
                log.warn("Detected an Intel MKL implementation of fftw resident in memory. Directly loading libfftw3.");
                dlhandle = dlmopen(LM_ID_NEWLM, dlInfo.dli_fname, RTLD_LOCAL | RTLD_DEEPBIND | RTLD_LAZY);
            }

            return dlhandle;
        #else
            return NULL;
        #endif
    }

    friend std::auto_ptr<libfftw3_t>;
};

#undef DECLARE_FFTW
#undef BIND_FFTW

#endif // !LSST_MEAS_BASE_libfftw_h_INCLUDED
