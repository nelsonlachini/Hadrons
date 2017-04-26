/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/Quark.hpp

Copyright (C) 2015
Copyright (C) 2016

Author: Antonin Portelli <antonin.portelli@me.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_Quark_hpp_
#define Hadrons_Quark_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 * 5D -> 4D and 4D -> 5D conversions.                                         *
 ******************************************************************************/
template<class vobj> // Note that 5D object is modified.
inline void make_4D(Lattice<vobj> &in_5d, Lattice<vobj> &out_4d, int Ls)
{
    axpby_ssp_pminus(in_5d, 0., in_5d, 1., in_5d, 0, 0);
    axpby_ssp_pplus(in_5d, 1., in_5d, 1., in_5d, 0, Ls-1);
    ExtractSlice(out_4d, in_5d, 0, 0);
}

template<class vobj>
inline void make_5D(Lattice<vobj> &in_4d, Lattice<vobj> &out_5d, int Ls)
{
    out_5d = zero;
    InsertSlice(in_4d, out_5d, 0, 0);
    InsertSlice(in_4d, out_5d, Ls-1, 0);
    axpby_ssp_pplus(out_5d, 0., out_5d, 1., out_5d, 0, 0);
    axpby_ssp_pminus(out_5d, 0., out_5d, 1., out_5d, Ls-1, Ls-1);
}

/******************************************************************************
 *                               TQuark                                       *
 ******************************************************************************/
class QuarkPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(QuarkPar,
                                    std::string, source,
                                    std::string, solver);
};

template <typename FImpl>
class TQuark: public Module<QuarkPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TQuark(const std::string name);
    // destructor
    virtual ~TQuark(void) = default;
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
    SolverFn     *solver_{nullptr};
};

MODULE_REGISTER(Quark, TQuark<FIMPL>);

/******************************************************************************
 *                          TQuark implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TQuark<FImpl>::TQuark(const std::string name)
: Module(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TQuark<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().solver};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TQuark<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_5d"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TQuark<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().solver);
    env().template registerLattice<PropagatorField>(getName());
    if (Ls_ > 1)
    {
        env().template registerLattice<PropagatorField>(getName() + "_5d", Ls_);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TQuark<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
                 << std::endl;
    
    FermionField    source(env().getGrid(Ls_)), sol(env().getGrid(Ls_)),
                    tmp(env().getGrid());
    std::string     propName = (Ls_ == 1) ? getName() : (getName() + "_5d");
    PropagatorField &prop    = *env().template createLattice<PropagatorField>(propName);
    PropagatorField &fullSrc = *env().template getObject<PropagatorField>(par().source);
    SolverFn        &solver  = *env().template getObject<SolverFn>(par().solver);
    if (Ls_ > 1)
    {
        env().template createLattice<PropagatorField>(getName());
    }
    
    LOG(Message) << "Inverting using solver '" << par().solver
                 << "' on source '" << par().source << "'" << std::endl;
    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < Nc; ++c)
    {
        LOG(Message) << "Inversion for spin= " << s << ", color= " << c
        << std::endl;
        // source conversion for 4D sources
        if (!env().isObject5d(par().source))
        {
            if (Ls_ == 1)
            {
                PropToFerm(source, fullSrc, s, c);
            }
            else
            {
                PropToFerm(tmp, fullSrc, s, c);
                make_5D(tmp, source, Ls_);
            }
        }
        // source conversion for 5D sources
        else
        {
            if (Ls_ != env().getObjectLs(par().source))
            {
                HADRON_ERROR("Ls mismatch between quark action and source");
            }
            else
            {
                PropToFerm(source, fullSrc, s, c);
            }
        }
        sol = zero;
        solver(sol, source);
        FermToProp(prop, sol, s, c);
        // create 4D propagators from 5D one if necessary
        if (Ls_ > 1)
        {
            PropagatorField &p4d =
                *env().template getObject<PropagatorField>(getName());
            make_4D(sol, tmp, Ls_);
            FermToProp(p4d, tmp, s, c);
        }
    }
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Quark_hpp_
