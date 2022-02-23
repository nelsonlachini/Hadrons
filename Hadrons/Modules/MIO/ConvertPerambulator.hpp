/*
 * ConvertPerambulator.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 *  Author: Felix Erben <ferben@ed.ac.uk>
 *  Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <dc-erbe1@tesseract-login1.ib0.sgi.cluster.dirac.ed.ac.uk>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: ferben <ferben@debian.felix.com>
 * Author: ferben <ferben@localhost.localdomain>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */

#ifndef Hadrons_MIO_ConvertPerambulator_hpp_
#define Hadrons_MIO_ConvertPerambulator_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MIO)

/******************************************************************************
 *    ConvertPerambulator (in principle, only works for exact distillation    *
 ******************************************************************************/

class ConvertPerambulatorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConvertPerambulatorPar,
                                        std::string, perambFileName,
                                        std::string, newPerambFileName,
                                        std::string, distilNoise,
                                        std::string, lapEigenPack,
                                        std::string, timeSources);
};

template <typename FImpl>
class TConvertPerambulator: public Module<ConvertPerambulatorPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TConvertPerambulator(const std::string name);
    // destructor
    virtual ~TConvertPerambulator(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ConvertPerambulator, TConvertPerambulator<FIMPL>, MIO);

/******************************************************************************
 *                 TConvertPerambulator implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TConvertPerambulator<FImpl>::TConvertPerambulator(const std::string name) : Module<ConvertPerambulatorPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TConvertPerambulator<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().lapEigenPack, par().distilNoise};
    return {in};
}

template <typename FImpl>
std::vector<std::string> TConvertPerambulator<FImpl>::getOutput(void)
{
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TConvertPerambulator<FImpl>::setup(void)
{
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().distilNoise);
    int nNoise = dilNoise.size();        
    int nVec = dilNoise.getNl();        
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        
    const int  Nt{env().getDim(Tdir)};
    int nSourceT;
    std::string sourceT = par().timeSources;
    nSourceT = MDistil::verifyTimeSourcesInput(sourceT,nDT);
    
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);


    envCreate(MDistil::PerambTensor, getName(), 1, Nt, nVec, nDL, nNoise, nSourceT, nDS);
    envTmp(MDistil::PerambIndexTensor, "PerambTmp", 1, Nt, nVec, nDL, nNoise, nDS);
    envTmp(std::vector<typename DistillationNoise<FImpl>::LapPack>,  "epack_3d", 1, Nt);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TConvertPerambulator<FImpl>::execute(void)
{
    auto &epack_4d = envGet(typename DistillationNoise<FImpl>::LapPack, par().lapEigenPack);
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const unsigned int Nt_first = gridHD->LocalStarts()[gridHD->Nd() - 1];
    const unsigned int Nt_local = gridHD->LocalDimensions()[gridHD->Nd() - 1];
    
    auto &perambulator = envGet(MDistil::PerambTensor, getName());
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().distilNoise);
    int nNoise = dilNoise.size();        
    int nVec = dilNoise.getNl();        
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        
    const int  Nt{env().getDim(Tdir)};

    std::string sourceT = par().timeSources;
    int nSourceT;
    std::vector<int> invT;
    nSourceT = MDistil::getSourceTimesFromInput(sourceT,nDT,dilNoise,invT);    
    perambulator.MetaData.timeSources = invT;
    
    envGetTmp(MDistil::PerambIndexTensor, PerambTmp);
    envGetTmp(std::vector<typename DistillationNoise<FImpl>::LapPack>, epack_3d);   // Eigenpack for each timeslice

    for (unsigned int t = Nt_first; t < Nt_first + Nt_local; t++)
    {
        epack_3d[t].resize(epack_4d.evec.size(),gridLD);
        for (int i=0;i<nVec;i++)
        {
            ExtractSliceLocal(epack_3d[t].evec[i],epack_4d.evec[i],0,t-Nt_first,Tdir); // switch to 3d object
        }
    }

    for (int dt = 0; dt < Nt; dt++)
    {
        std::vector<int>::iterator it = std::find(std::begin(invT), std::end(invT), dt);
        //skip dilution indices which are not in invT
        if(it == std::end(invT))
        {
            continue;
        }
        LOG(Message) <<  "reading perambulator dt= " << dt << std::endl;
        int idt = it - std::begin(invT);
        std::string sPerambName {par().perambFileName};
        sPerambName.append("/iDT_");
        sPerambName.append(std::to_string(dt));
        sPerambName.append(".");
        sPerambName.append(std::to_string(vm().getTrajectory()));
        if(gridHD->IsBoss())
            PerambTmp.read(sPerambName.c_str());

        std::string sNewPerambName {par().newPerambFileName};
        Hadrons::mkdir(sNewPerambName);
        sNewPerambName.append("/iDT_");
        sNewPerambName.append(std::to_string(dt));
        sNewPerambName.append(".");
        sNewPerambName.append(std::to_string(vm().getTrajectory()));

        // multiplying perambulator by evec phases
        for (int t = 0; t < Nt; t++)
        {
            auto grid = epack_3d[t].evec[0].Grid();
            Coordinate siteFirst(grid->Nd(),0);

            ColourVector cv0_old_t;
            peekSite(cv0_old_t, epack_3d[t].evec[0], siteFirst);
            const std::complex<Real> cplx0_old_t{cv0_old_t()()(0).real(), cv0_old_t()()(0).imag()};
            if( cplx0_old_t.imag() == 0 )
                ;
            else
            {
                const Real cplx0_mag_old_t{ std::abs(cplx0_old_t) };
                const std::complex<Real> std_phase_old_t{std::conj(cplx0_old_t/cplx0_mag_old_t)};
                const Grid::Complex phase_old_t{std_phase_old_t.real(),std_phase_old_t.imag()};

                ColourVector cv0_old_dt;
                peekSite(cv0_old_dt, epack_3d[dt].evec[0], siteFirst);
                const std::complex<Real> cplx0_old_dt{cv0_old_dt()()(0).real(), cv0_old_dt()()(0).imag()};
                if( cplx0_old_t.imag() == 0 )
                    ;
                else
                {
                    const Real cplx0_mag_old_dt{ std::abs(cplx0_old_dt) };
                    const std::complex<Real> std_phase_old_dt{std::conj(cplx0_old_dt/cplx0_mag_old_dt)};
                    const Grid::Complex phase_old_dt{std_phase_old_dt.real(),-std_phase_old_dt.imag()};

                    const Grid::Complex phi_old_conj  = phase_old_t * phase_old_dt;

                    for( int dk = 0 ; dk < nVec ; dk++ )
                    {
                        ColourVector cv0_dk;
                        peekSite(cv0_dk, epack_3d[dt].evec[dk], siteFirst);
                        const std::complex<Real> cplx0_dk{cv0_dk()()(0).real(), cv0_dk()()(0).imag()};
                        if( cplx0_dk.imag() == 0 )
                            ;
                        else
                        {
                            const Real cplx0_mag_dk{ std::abs(cplx0_dk) };
                            const std::complex<Real> std_phase_dk{std::conj(cplx0_dk/cplx0_mag_dk)};
                            const Grid::Complex phase_dk{std_phase_dk.real(),std_phase_dk.imag()};

                            for (int ivec = 0; ivec < nVec; ivec++)
                            {
                                ColourVector cv0_ivec;
                                peekSite(cv0_ivec, epack_3d[t].evec[ivec], siteFirst);
                                const std::complex<Real> cplx0_ivec{cv0_ivec()()(0).real(), cv0_ivec()()(0).imag()};
                                if( cplx0_ivec.imag() == 0 )
                                    ;
                                else
                                {
                                    const Real cplx0_mag_ivec{ std::abs(cplx0_ivec) };
                                    const std::complex<Real> std_phase_ivec{std::conj(cplx0_ivec/cplx0_mag_ivec)};
                                    const Grid::Complex phase_ivec{std_phase_ivec.real(),-std_phase_ivec.imag()};

                                    const Grid::Complex phi_new  = phase_dk * phase_ivec;

                                    for (int in = 0; in < nNoise; in++)
                                    for (int ids = 0; ids < nDS; ids++)
                                    {
                                        PerambTmp.tensor(t,ivec,dk,in,ids) = PerambTmp.tensor(t,ivec,dk,in,ids) * phi_new * phi_old_conj;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if(gridHD->IsBoss())
            PerambTmp.write(sNewPerambName.c_str());
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
