#ifndef Hadrons_MIO_LoadPerambulatorTimeSourceMetadata_hpp_
#define Hadrons_MIO_LoadPerambulatorTimeSourceMetadata_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/NamedTensorTimeSourceMetadata.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LoadPerambulatorTimeSourceMetadata                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadPerambulatorTimeSourceMetadataPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadPerambulatorTimeSourceMetadataPar,
                                        std::string, oldPerambFileName,
                                        std::string, perambFileName,
                                        int, ti,
                                        int, li,
                                        int, nVec,
                                        std::string, timeSources,
                                        std::string, nNoise);
};

template <typename FImpl>
class TLoadPerambulatorTimeSourceMetadata: public Module<LoadPerambulatorTimeSourceMetadataPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TLoadPerambulatorTimeSourceMetadata(const std::string name);
    // destructor
    virtual ~TLoadPerambulatorTimeSourceMetadata(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    int nnoise_=1;
    std::vector<int> ts_;
    int ntInv_=1;
};

MODULE_REGISTER_TMP(LoadPerambulatorTimeSourceMetadata, TLoadPerambulatorTimeSourceMetadata<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadPerambulatorTimeSourceMetadata implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadPerambulatorTimeSourceMetadata<FImpl>::TLoadPerambulatorTimeSourceMetadata(const std::string name)
: Module<LoadPerambulatorTimeSourceMetadataPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadPerambulatorTimeSourceMetadata<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadPerambulatorTimeSourceMetadata<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPerambulatorTimeSourceMetadata<FImpl>::setup(void)
{
    if(par().nNoise!="")
    {
        nnoise_ = stoi(par().nNoise);
    }
    ts_ = strToVec<int>(par().timeSources);
    ntInv_ = ts_.size();
    
    // std::cout << "Using sizes " << env().getDim(Tdir) << " " << par().nVec << " " << par().li << " " << nnoise_ << " " << ntInv_ << " " << Ns << std::endl;
    envCreate(MDistil::PerambTensorTimeSourceMetadata, getName()+"TimeSourceMetadata", 1, env().getDim(Tdir),
                     par().nVec, par().li, nnoise_, ntInv_, Ns);
    envTmp(MDistil::PerambIndexTensor, "PerambDT", 1, env().getDim(Tdir), par().nVec, par().li, nnoise_, Ns);
    // std::cin.get();
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPerambulatorTimeSourceMetadata<FImpl>::execute(void)
{
    GridCartesian * grid4d = envGetGrid(FermionField);
    auto &peramb_old = envGet(MDistil::PerambTensorTimeSourceMetadata, getName()+"TimeSourceMetadata");
    std::string sOldPerambName{ par().oldPerambFileName };
    sOldPerambName.append( 1, '.' );
    sOldPerambName.append( std::to_string( vm().getTrajectory() ) );
    peramb_old.read(sOldPerambName.c_str());

    envGetTmp(MDistil::PerambIndexTensor, PerambDT);
    int idt=0;
    for(auto dt : ts_)
    {
        std::string sNewPerambName {par().perambFileName};
        sNewPerambName.append("/iDT_");
        sNewPerambName.append(std::to_string(dt));
        sNewPerambName.append(".");
        sNewPerambName.append(std::to_string(vm().getTrajectory()));
        makeFileDir(sNewPerambName, grid4d);
        for(int t=0         ; t<env().getDim(Tdir)  ; t++)
        for(int ivec=0      ; ivec<par().nVec       ; ivec++)
        for(int dl=0        ; dl<par().li           ; dl++)
        for(int inoise=0    ; inoise<nnoise_        ; inoise++)
        for(int ds=0        ; ds<Ns                 ; ds++)
        {
            // std::cout << t << " " << ivec << " " << dl << " " << inoise << " " << idt << " " << ds << " " << peramb_old.tensor(t,ivec,dl,inoise,idt,ds) << std::endl;
            PerambDT.tensor(t,ivec,dl,inoise,ds) = peramb_old.tensor(t,ivec,dl,inoise,idt,ds);
        }
        PerambDT.MetaData.timeDilutionIndex = dt;
        PerambDT.write(sNewPerambName.c_str());
        // std::cin.get();
        idt++;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadPerambulatorTimeSourceMetadata_hpp_
