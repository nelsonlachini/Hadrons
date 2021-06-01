#ifndef Hadrons_MIO_LoadPerambulatorNoMetadata_hpp_
#define Hadrons_MIO_LoadPerambulatorNoMetadata_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/NamedTensorNoMetadata.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LoadPerambulatorNoMetadata                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadPerambulatorNoMetadataPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadPerambulatorNoMetadataPar,
                                        std::string, oldPerambFileName,
                                        std::string, perambFileName,
                                        int, ti,
                                        int, li,
                                        int, nVec,
                                        std::string, timeSources,
                                        std::string, nNoise);
};

template <typename FImpl>
class TLoadPerambulatorNoMetadata: public Module<LoadPerambulatorNoMetadataPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TLoadPerambulatorNoMetadata(const std::string name);
    // destructor
    virtual ~TLoadPerambulatorNoMetadata(void) {};
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

MODULE_REGISTER_TMP(LoadPerambulatorNoMetadata, TLoadPerambulatorNoMetadata<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadPerambulatorNoMetadata implementation                 *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadPerambulatorNoMetadata<FImpl>::TLoadPerambulatorNoMetadata(const std::string name)
: Module<LoadPerambulatorNoMetadataPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadPerambulatorNoMetadata<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadPerambulatorNoMetadata<FImpl>::getOutput(void)
{
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPerambulatorNoMetadata<FImpl>::setup(void)
{
    if(par().nNoise!="")
    {
        nnoise_ = stoi(par().nNoise);
    }
    ts_ = strToVec<int>(par().timeSources);
    ntInv_ = ts_.size();
    
    // std::cout << "Using sizes " << env().getDim(Tdir) << " " << par().nVec << " " << par().li << " " << nnoise_ << " " << ntInv_ << " " << Ns << std::endl;
    envCreate(MDistil::PerambTensorNoMetadata, getName()+"NoMetadata", 1, env().getDim(Tdir),
                     par().nVec, par().li, nnoise_, ntInv_, Ns);
    envTmp(MDistil::PerambIndexTensor, "PerambDT", 1, env().getDim(Tdir), par().nVec, par().li, nnoise_, Ns);
    // std::cin.get();
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadPerambulatorNoMetadata<FImpl>::execute(void)
{
    GridCartesian * grid4d = envGetGrid(FermionField);
    auto &peramb_old = envGet(MDistil::PerambTensorNoMetadata, getName()+"NoMetadata");
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

#endif // Hadrons_MIO_LoadPerambulatorNoMetadata_hpp_
