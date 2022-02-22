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

    envCreate(MDistil::PerambTensor, getName(), 1, Nt, nVec, nDL, nNoise, nSourceT, nDS);
    envTmp(MDistil::PerambIndexTensor, "PerambTmp", 1, Nt, nVec, nDL, nNoise, nDS);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TConvertPerambulator<FImpl>::execute(void)
{
    auto &epack = envGet(typename DistillationNoise<FImpl>::LapPack, par().lapEigenPack);
    GridCartesian *grid            = envGetGrid(ColourVectorField);
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
        PerambTmp.read(sPerambName.c_str());

        std::string sNewPerambName {par().newPerambFileName};
        Hadrons::mkdir(sNewPerambName);
        sNewPerambName.append("/iDT_");
        sNewPerambName.append(std::to_string(dt));
        sNewPerambName.append(".");
        sNewPerambName.append(std::to_string(vm().getTrajectory()));

        // perambulator.tensor(t,ivec,idl,in,idt,ids) = PerambTmp.tensor(t,ivec,idl,in,ids);

        // multiplying perambulator by evec phases
        auto evec = epack.evec;
        Coordinate siteFirst(grid->Nd(),0);
    
        for( int dk = 0 ; dk < evec.size() ; dk++ )
        {
            ColourVector cv0_dk;
            peekSite(cv0_dk, evec[dk], siteFirst);
            const std::complex<Real> cplx0_dk{cv0_dk()()(0).real(), cv0_dk()()(0).imag()};
            if( cplx0_dk.imag() == 0 )
                // LOG(Message) << "RotateEigen() : Vector " << dk <<  " Site 0 : " << cplx0_dk << " => already meets phase convention" << std::endl;
                ;
            else
            {
                const Real cplx0_mag_dk{ std::abs(cplx0_dk) };
                const std::complex<Real> std_phase_dk{std::conj(cplx0_dk/cplx0_mag_dk)};
                // LOG(Message) << "RotateEigen() : Vector " << dk <<  " Site 0 : |" << cplx0_dk << "|=" << cplx0_mag_dk
                        // << " => phase=" << (std::arg(std_phase_dk) / M_PI) << " pi" << std::endl;
                const Grid::Complex phase_dk{std_phase_dk.real(),std_phase_dk.imag()};

                for (int ivec = 0; ivec < nVec; ivec++)
                {
                    ColourVector cv0_ivec;
                    peekSite(cv0_ivec, evec[ivec], siteFirst);
                    const std::complex<Real> cplx0_ivec{cv0_ivec()()(0).real(), cv0_ivec()()(0).imag()};
                    if( cplx0_ivec.imag() == 0 )
                        ;
                    else
                    {
                        const Real cplx0_mag_ivec{ std::abs(cplx0_ivec) };
                        const std::complex<Real> std_phase_ivec{std::conj(cplx0_ivec/cplx0_mag_ivec)};
                        const Grid::Complex phase_ivec{std_phase_ivec.real(),-std_phase_ivec.imag()};

                        for (int t = 0; t < Nt; t++)
                        for (int in = 0; in < nNoise; in++)
                        for (int ids = 0; ids < nDS; ids++)
                        {
                            PerambTmp.tensor(t,ivec,dk,in,ids) = PerambTmp.tensor(t,ivec,dk,in,ids) * phase_dk * phase_ivec;
                        }
                    }
                }
            }
        }
        PerambTmp.write(sNewPerambName.c_str());
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
