/*
 *	This file is part of LCQPow.
 *
 *	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2021 by Jonas Hall et al.
 *
 *	LCQPow is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	LCQPow is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with LCQPow; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "SubsolverOSQP.hpp"

extern "C" {
    #include <osqp.h>
}


namespace LCQPow {
    SubsolverOSQP::SubsolverOSQP( ) {
        data = NULL;
        settings = NULL;
        work = NULL;
     }


    SubsolverOSQP::SubsolverOSQP(   const csc* const _H, const csc* const _A)
    {
        // Store dimensions
        nV = _H->n;
        nC = _A->m;

        // Allocate memory for settings
        settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

        // Copy matrices
        H = Utilities::copyCSC(_H, true);
        A = Utilities::copyCSC(_A);

        // Define solver settings
        osqp_set_default_settings(settings);
        settings->eps_prim_inf = Utilities::ZERO;
        settings->verbose = false;
        settings->polish = true;
    }


    SubsolverOSQP::SubsolverOSQP(const SubsolverOSQP& rhs)
    {
        copy( rhs );
    }


    SubsolverOSQP::~SubsolverOSQP()
    {
        clear();
    }


    void SubsolverOSQP::clear()
    {
        osqp_cleanup(work);

        if (data != 0) {
            if (data->l != 0) {
                free(data->l);
                data->l = NULL;
            }

            if (data->u != 0) {
                free(data->u);
                data->u = NULL;
            }

            if (data->q != 0) {
                free(data->q);
                data->q = NULL;
            }

            c_free(data);
            data = NULL;
        }

        if (H != 0) {
            Utilities::ClearSparseMat(H);
            H = NULL;
        }

        if (A != 0) {
            Utilities::ClearSparseMat(A);
            A = NULL;
        }

        if (settings != 0) {
            c_free(settings);
            settings = NULL;
        }
    }


    SubsolverOSQP& SubsolverOSQP::operator=(const SubsolverOSQP& rhs)
    {
        if (this != &rhs) {
            copy( rhs );
        }

        return *this;
    }


    void SubsolverOSQP::setOptions( OSQPSettings* _settings )
    {
        settings = copy_settings(_settings);
    }


    void SubsolverOSQP::setPrintlevl( bool verbose )
    {
        if (settings != 0)
            settings->verbose = verbose;
    }


    ReturnValue SubsolverOSQP::solve(   bool initialSolve, int& iterations,
                                        const double* const _g,
                                        const double* const _lbA, const double* const _ubA,
                                        const double* const x0, const double* const y0,
                                        const double* const _lb, const double* const _ub )
    {
        // Make sure that lb and ub are null pointers, as OSQP does not handle box constraints
        if (_lb != 0 || _ub != 0) {
            return ReturnValue::INVALID_OSQP_BOX_CONSTRAINTS;
        }

        // Setup workspace on initial solve
        if (initialSolve) {

            double* l = (double*)malloc((size_t)nC*sizeof(double));
            double* u = (double*)malloc((size_t)nC*sizeof(double));
            double* g = (double*)malloc((size_t)nV*sizeof(double));
            memcpy(l, _lbA, (size_t)nC*sizeof(double));
            memcpy(u, _ubA, (size_t)nC*sizeof(double));
            memcpy(g, _g, (size_t)nV*sizeof(c_float));

            data = (OSQPData *)c_malloc(sizeof(OSQPData));
            data->n = nV;
            data->m = nC;
            data->P = H;
            data->A = A;
            data->q = g;
            data->l = l;
            data->u = u;
            osqp_setup(&work, data, settings);
        } else {
            // Update linear cost and bounds
            osqp_update_lin_cost(work, _g);
            osqp_update_bounds(work, _lbA, _ubA);
        }

        if (work == 0) {
            return ReturnValue::OSQP_WORKSPACE_NOT_SET_UP;
        }

        // Solve Problem
        int exitflag = osqp_solve(work);

        // Get number of iterations
        iterations = work->info->iter;

        // Either pass error
        if (exitflag != 0)
            return ReturnValue::SUBPROBLEM_SOLVER_ERROR;

        // Or pass successful return
        return ReturnValue::SUCCESSFUL_RETURN;
    }


    void SubsolverOSQP::getSolution( double* x, double* y )
    {
        OSQPSolution *sol(work->solution);

        if (sol->x != 0) {
            memcpy(x, sol->x, (size_t)nV*(sizeof(double)));
        }

        // Copy duals with negative sign
        for (int i = 0; i < nC; i++) {
            y[i] = -sol->y[i];
        }
    }


    void SubsolverOSQP::copy(const SubsolverOSQP& rhs)
    {
        clear();

        nV = rhs.nV;
        nC = rhs.nC;

        H = copy_csc_mat(rhs.H);
        A = copy_csc_mat(rhs.A);

        if (rhs.settings != 0) {
            settings = copy_settings(rhs.settings);
        }

        if (rhs.data != 0) {
            double* l = (double*)malloc((size_t)nC*sizeof(double));
            double* u = (double*)malloc((size_t)nC*sizeof(double));
            double* g = (double*)malloc((size_t)nV*sizeof(double));
            memcpy(l, rhs.data->l, (size_t)nC*sizeof(double));
            memcpy(u, rhs.data->u, (size_t)nC*sizeof(double));
            memcpy(g, rhs.data->q, (size_t)nV*sizeof(c_float));

            data = (OSQPData *)c_malloc(sizeof(OSQPData));
            data->n = nV;
            data->m = nC;
            data->P = H;
            data->A = A;
            data->q = g;
            data->l = l;
            data->u = u;

            osqp_setup(&work, data, settings);
        }
    }
}