#ifndef Hadrons_MDistil_DistilMesonField_hpp_
#define Hadrons_MDistil_DistilMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MDistil/Distil.hpp>
#include <Hadrons/A2AMatrix.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DistilMesonField                                 *
 * Eliminates DistilVectors module. Receives LapH eigenvectors and 
 * perambulator/noise (as left/right fields). Computes MesonFields by 
 * block(and chunking it) and save them to H5 file.
 * 
 * For now, do not load anything from disk. Trying phi-phi case 
 * with full-dilution.
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MDistil)

class DistilMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra, gamma,
                                    std::vector<int>,   noise_pair,
                                    std::vector<std::vector<int>>, dilution_left,
                                    std::vector<std::vector<int>>, dilution_right);
};

class DistilMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldPar,
                                    std::string,    OutputStem,
                                    std::string,    MesonFieldCase,
                                    std::string,    LapEvec,
                                    std::string,    LeftInput,
                                    std::string,    RightInput,
                                    std::string,    LeftDPar,
                                    std::string,    RightDPar,
                                    std::vector<std::string>, NoisePairs,
                                    std::vector<std::string>, SourceTimesLeft,
                                    std::vector<std::string>, SourceTimesRight,
                                    std::string, Gamma,
                                    std::vector<std::string>, Momenta,
                                    int, BlockSize,
                                    int, CacheSize)
};

template <typename FImpl>
class TDistilMesonField: public Module<DistilMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TDistilMesonField(const std::string name);
    // destructor
    virtual ~TDistilMesonField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
protected:
    std::unique_ptr<GridCartesian> grid3d;
private:
    std::map<std::string,std::string>   dmf_case;
    // std::string                        momphName_;
    std::vector<Gamma::Algebra>         gamma_;
    std::vector<std::vector<RealF>>     momenta_;
    int                                 nExt_;
    int                                 nStr_;
    Vector<ComplexF>                    buf_;
    Vector<Complex>                     cache_;
    int                                 nt_nonzero_;
    std::vector<std::vector<int>>       noisePairs_;           // read from extermal object (diluted noise class)
    std::map<std::string, std::string>                      lrInput_;
    std::map<std::string, std::vector<std::vector<int>>>    lrSourceTimes_;
    std::string outputMFStem;
    bool                               hasPhase_{false};
};

MODULE_REGISTER_TMP(DistilMesonField, TDistilMesonField<FIMPL>, MDistil);

/******************************************************************************
 *                 TDistilMesonField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilMesonField<FImpl>::TDistilMesonField(const std::string name)
: Module<DistilMesonFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilMesonField<FImpl>::getInput(void)
{   
    return{par().LapEvec, par().LeftInput, par().RightInput, par().LeftDPar, par().RightDPar};
}

template <typename FImpl>
std::vector<std::string> TDistilMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}


// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonField<FImpl>::setup(void)
{
    outputMFStem = par().OutputStem;

    if(par().MesonFieldCase=="phi-phi" || par().MesonFieldCase=="phi-rho" || par().MesonFieldCase=="rho-phi" || par().MesonFieldCase=="rho-rho")
    {
        dmf_case.emplace("left",par().MesonFieldCase.substr(0,3));   //left
        dmf_case.emplace("right",par().MesonFieldCase.substr(4,3));  //right
    }
    else
    {
        HADRONS_ERROR(Argument,"Bad meson field case");
    }
    LOG(Message) << "Meson field case " << dmf_case.at("left") << "-" << dmf_case.at("right") << std::endl;
    const std::vector<std::string> lrPair = {"left","right"};

    // parse source times lr
    std::vector<std::vector<int>>       stL;      // should this come from perambulator in the future? no only if we want to use subsets of the inversions we have
    stL.clear();
    for(auto &stl : par().SourceTimesLeft)
        stL.push_back(strToVec<int>(stl));
        
    std::vector<std::vector<int>>       stR;
    stR.clear();
    for(auto &str : par().SourceTimesRight)
        stR.push_back(strToVec<int>(str));
    // outermost dimension is the time-dilution index, innermost one are the non-zero source timeslices
    // in phi-phi save all timeslices, but in the other cases save only the non-zero  ones...
    
    lrInput_       = {{"left",par().LeftInput},{"right",par().RightInput}};
    lrSourceTimes_ = {{"left",stL},{"right",stR}};  

    // parse and validate noise pairs
    noisePairs_.clear();
    for(auto &npair : par().NoisePairs)
    {
        noisePairs_.push_back(strToVec<int>(npair));
        std::map<std::string, int>    noiseMapTemp = { {"left", noisePairs_.back()[0]} , {"right",noisePairs_.back()[1]} };
        for(auto &side : lrPair){
            if(dmf_case.at(side)=="phi")
            {
                // turn this into macro?
                auto &inTensor = envGet(PerambTensor , lrInput_.at(side));
                if( noiseMapTemp.at(side) > inTensor.tensor.dimensions().at(3) )
                {
                    HADRONS_ERROR(Size,"Noise pair element " + std::to_string(noiseMapTemp.at(side)) + ">" +std::to_string(inTensor.tensor.dimensions().at(3)) + " unavailable on input tensor");
                }
            }
            else
            {
                auto &inTensor = envGet(NoiseTensor , lrInput_.at(side));
                if( noiseMapTemp.at(side) > inTensor.tensor.dimensions().at(3) )
                {
                    HADRONS_ERROR(Size,"Noise pair element " + std::to_string(noiseMapTemp.at(side)) + ">" +std::to_string(inTensor.tensor.dimensions().at(3)) + " unavailable on input tensor");
                }
            }
        }
    }

    // momenta and gamma parse
    momenta_.clear();
    for(auto &p_string : par().Momenta)
    {
        auto p = strToVec<RealF>(p_string);

        if (p.size() != env().getNd() - 1)
        {
            HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size())
                                + " components instead of " 
                                + std::to_string(env().getNd() - 1));
        }
        momenta_.push_back(p);
    }
    if (par().Gamma == "all")
    {
        gamma_ = {
            Gamma::Algebra::Gamma5,
            Gamma::Algebra::Identity,    
            Gamma::Algebra::GammaX,
            Gamma::Algebra::GammaY,
            Gamma::Algebra::GammaZ,
            Gamma::Algebra::GammaT,
            Gamma::Algebra::GammaXGamma5,
            Gamma::Algebra::GammaYGamma5,
            Gamma::Algebra::GammaZGamma5,
            Gamma::Algebra::GammaTGamma5,
            Gamma::Algebra::SigmaXY,
            Gamma::Algebra::SigmaXZ,
            Gamma::Algebra::SigmaXT,
            Gamma::Algebra::SigmaYZ,
            Gamma::Algebra::SigmaYT,
            Gamma::Algebra::SigmaZT
        };
    }
    else
    {
        gamma_ = strToVec<Gamma::Algebra>(par().Gamma);
    }

    // hard-coded dilution scheme (assuming dpL=dpR for now)
    const DistilParameters &dpL = envGet(DistilParameters, par().LeftDPar);
    const DistilParameters &dpR = envGet(DistilParameters, par().RightDPar);

    //matrix sets
    nExt_ = momenta_.size(); // * nNoiseLeft * nNoiseRight;   // remember to take into account the noise pair in the filename later
    nStr_ = gamma_.size();

    // compute nt_nonzero_ (<nt), the number of non-zero timeslices when there's at least one rho involved
    nt_nonzero_ = 0;
    if(par().MesonFieldCase=="rho-rho" || par().MesonFieldCase=="rho-phi")
        for(auto elem : stL)
            elem.size() > nt_nonzero_ ? nt_nonzero_ = elem.size() : 0;      //get the highest possible nt_nonzero_ from stL
    else if(par().MesonFieldCase=="phi-rho")
        for(auto elem : stR)
            elem.size() > nt_nonzero_ ? nt_nonzero_ = elem.size() : 0;
    else
        nt_nonzero_ = env().getDim(Tdir);
    
    // std::cout << "Smallest time extension = " << nt_nonzero_ << std::endl << std::endl;
    buf_.resize(nExt_*nStr_*nt_nonzero_*par().BlockSize*par().BlockSize);
    cache_.resize(nExt_*nStr_*env().getDim(Tdir)*par().CacheSize*par().CacheSize);
    
    // 3d grid (as a 4d one with collapsed time dimension)
    const unsigned int nd{env().getNd()};
    GridCartesian * const grid4d{env().getGrid()};
    Coordinate latt_size_3          = grid4d->_gdimensions;
    latt_size_3[nd-1]               = 1;
    Coordinate simd_layout_3        = GridDefaultSimd(nd-1, vComplex::Nsimd());
    simd_layout_3.push_back( 1 );
    Coordinate mpi_layout_3         = grid4d->_processors;
    mpi_layout_3[nd-1]              = 1;
    grid3d.reset( new GridCartesian(latt_size_3, simd_layout_3, mpi_layout_3, *grid4d) );

    envTmp(FermionField,                    "fermion3dtmp",         1, grid3d.get());
    envTmp(ColourVectorField,               "fermion3dtmp_nospin",  1, grid3d.get());
    envTmp(ColourVectorField,               "evec3d",               1, grid3d.get());
    envTmp(std::vector<FermionField>,       "left",                 1, stL.size()*dpL.LI*dpL.SI,   envGetGrid(FermionField));
    envTmp(std::vector<FermionField>,       "right",                1, stL.size()*dpR.LI*dpR.SI,   envGetGrid(FermionField));
    envTmpLat(ComplexField, "coor");
    envCache(std::vector<ComplexField>,     "phasename",            1, momenta_.size(), envGetGrid(ComplexField));
    envTmpLat(FermionField,                 "fermion4dtmp");

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilMesonField<FImpl>::execute(void)
{
    // temps
    envGetTmp(FermionField,                 fermion3dtmp);
    envGetTmp(ColourVectorField,            fermion3dtmp_nospin);
    envGetTmp(FermionField,                 fermion4dtmp);
    envGetTmp(ColourVectorField,            evec3d);
    envGetTmp(std::vector<FermionField>,    left);
    envGetTmp(std::vector<FermionField>,    right);

    int blockSize_{par().BlockSize};
    int cacheSize_{par().CacheSize};
    const int nt{env().getDim(Tdir)};
    const unsigned int nd{env().getNd()};
    const int Ntlocal{ env().getGrid()->LocalDimensions()[Tdir] };
    const int Ntfirst{ env().getGrid()->LocalStarts()[Tdir] };
    int vol = 1;
    for(int v; v<Nd ; v++)
    {
        vol *= env().getGrid()->GlobalDimensions()[v];
    }
    
    auto &epack = envGet(LapEvecs, par().LapEvec);
    auto &phase = envGet(std::vector<ComplexField>, "phasename");
    if (!hasPhase_)
    {
        startTimer("momentum phases");
        for (unsigned int j = 0; j < momenta_.size(); ++j)
        {
            Complex           i(0.0,1.0);
            std::vector<Real> p;

            envGetTmp(ComplexField, coor);
            phase[j] = Zero();
            for(unsigned int mu = 0; mu < momenta_[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                phase[j] = phase[j] + (momenta_[j][mu]/env().getDim(mu))*coor;
            }
            phase[j] = exp((Real)(2*M_PI)*i*phase[j]);
        }
        hasPhase_ = true;
        stopTimer("momentum phases");
    }

    // hard-coded dilution schem &  assuming dilution_left == dilution_right; other cases?...
    // replace by noise class
    const DistilParameters &dpL = envGet(DistilParameters, par().LeftDPar);
    const DistilParameters &dpR = envGet(DistilParameters, par().RightDPar);
    
    int nInversions             = MIN(dpL.inversions,dpR.inversions);   // should nInversions be greater or equal to SourceTimesLeft.size() == SourceTimesRight.size() always? I guess so
    assert(nInversions >= lrSourceTimes_.at("left").size());
    assert(nInversions >= lrSourceTimes_.at("right").size());

    //clean this out
    int LI = dpL.LI , SI = dpL.SI;
    int dilutionSize_LS = LI * SI;

    int nVec        = dpL.nvec;
    int nNoiseLeft  = dpL.nnoise;
    int nNoiseRight = dpR.nnoise;

    // do not use operator []!! similar but better way to do that? maybe map of pointers?
    std::map<std::string, std::vector<FermionField>&>       lrDistVector  = {{"left",left},{"right",right}};
    const std::vector<std::string> lrPair = {"left","right"};

    for(auto &inoise : noisePairs_)
    {
        //set up io object and metadata for all gamma/momenta
        std::vector<A2AMatrixIo<ComplexF>> matrixIoTable;
        DistilMesonFieldMetadata md;
        for(int iExt=0; iExt<nExt_; iExt++)
        for(int iStr=0; iStr<nStr_; iStr++)
        {
            // metadata;
            md.momentum = momenta_[iExt];
            md.gamma = gamma_[iStr];
            md.noise_pair = inoise;
            md.dilution_left = lrSourceTimes_.at("left");
            md.dilution_right = lrSourceTimes_.at("right");

            std::stringstream ss;
            ss << md.gamma << "_";
            for (unsigned int mu = 0; mu < md.momentum.size(); ++mu)
                ss << md.momentum[mu] << ((mu == md.momentum.size() - 1) ? "" : "_");
            std::string groupName = ss.str();

            //init file here (do not create dataset yet)
            //IO configuration for fixed test gamma and momentum
            
            std::string outputStem = outputMFStem + "/noise" + std::to_string(inoise[0]) + "_" + std::to_string(inoise[1]) + "/";
            Hadrons::mkdir(outputStem);
            A2AMatrixIo<ComplexF> matrixIo(outputStem+groupName+"_"+par().MesonFieldCase+".h5", groupName, nt_nonzero_, dilutionSize_LS, dilutionSize_LS);  // automatise name choice according to momenta_ and gamma_
            matrixIoTable.push_back(matrixIo);

            //initialize file with no outputName group (containing atributes of momentum and gamma) but no dataset inside
            if(env().getGrid()->IsBoss())
            {
                matrixIoTable.back().initFile(md);
            }
        }

        LOG(Message) << "Gamma:" << std::endl;
        LOG(Message) << gamma_ << std::endl;
        LOG(Message) << "Momenta:" << std::endl;
        LOG(Message) << momenta_ << std::endl;

        std::map<std::string,int&> lfNoise = {{"left",inoise[0]},{"right",inoise[1]}}; //do not use []

        // computation, still ignoring gamma5 hermiticity
        for(auto &side : lrPair)
        {
            //each dtL, dtR pair corresponds to a different dataset
            for (int dt = 0; dt<lrSourceTimes_.at(side).size(); dt++)         //loop over left time dilution index
            {
                // left computation
                for(int iD=0 ; iD<dilutionSize_LS ; iD++)
                {
                    int dk = iD%LI;
                    int ds = (iD - dk)/LI;
                    lrDistVector.at(side)[dt*dilutionSize_LS + iD] = Zero();
                    if(dmf_case.at(side)=="phi")
                    {
                        auto &inTensor = envGet(PerambTensor , lrInput_.at(side));
                        for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++)   //loop over (local) timeslices
                        {
                            fermion3dtmp = Zero();
                            for (int k = 0; k < nVec; k++)
                            {
                                ExtractSliceLocal(evec3d,epack.evec[k],0,t-Ntfirst,Tdir);
                                fermion3dtmp += evec3d * inTensor.tensor(t, k, dk, lfNoise.at(side), dt, ds);
                            }
                            InsertSliceLocal(fermion3dtmp,lrDistVector.at(side)[dt*dilutionSize_LS + iD],0,t-Ntfirst,Tdir);
                        }
                    }
                    else if(dmf_case.at(side)=="rho"){
                        auto &inTensor = envGet(NoiseTensor , lrInput_.at(side));
                        for(int t : lrSourceTimes_.at(side)[dt])
                        {
                            if (t >= Ntfirst && t < Ntfirst + Ntlocal)
                            {
                                fermion3dtmp = Zero();
                                for (int k = dk; k < nVec; k += LI)
                                {
                                    for (int s = ds; s < Ns; s += SI)
                                    {
                                        ExtractSliceLocal(evec3d,epack.evec[k],0,t-Ntfirst,Tdir);
                                        fermion3dtmp_nospin = evec3d * inTensor.tensor(lfNoise.at(side), t, k , s);
                                        fermion3dtmp = Zero();
                                        pokeSpin(fermion3dtmp,fermion3dtmp_nospin,s);
                                        fermion4dtmp = Zero();
                                        InsertSliceLocal(fermion3dtmp,fermion4dtmp,0,t-Ntfirst,Tdir);
                                        lrDistVector.at(side)[dt*dilutionSize_LS + iD] += fermion4dtmp;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // computing mesonfield blocks and saving to disk
        for (int dtL = 0; dtL < lrSourceTimes_.at("left").size() ; dtL++)
        for (int dtR = 0; dtR < lrSourceTimes_.at("right").size() ; dtR++)
        {
            if(!(par().MesonFieldCase=="rho-rho" && dtL!=dtR))
            {
                std::string datasetName = "dtL"+std::to_string(dtL)+"_dtR"+std::to_string(dtR);
                LOG(Message) << "Computing dilution block " << datasetName << std::endl;

                // int nblocki = left.size()/blockSize_ + (((left.size() % blockSize_) != 0) ? 1 : 0);
                // int nblockj = right.size()/blockSize_ + (((right.size() % blockSize_) != 0) ? 1 : 0);
                int nblocki = dilutionSize_LS/blockSize_ + (((dilutionSize_LS % blockSize_) != 0) ? 1 : 0);
                int nblockj = dilutionSize_LS/blockSize_ + (((dilutionSize_LS % blockSize_) != 0) ? 1 : 0);

                // loop over blocks in the current time-dilution block
                for(int iblock=0 ; iblock<dilutionSize_LS ; iblock+=blockSize_) //set according to memory size
                for(int jblock=0 ; jblock<dilutionSize_LS ; jblock+=blockSize_)
                {
                    int iblockSize = MIN(dilutionSize_LS-iblock,blockSize_);    // iblockSize is the size of the current block (indexed by i); N_i-i is the size of the eventual remainder block
                    int jblockSize = MIN(dilutionSize_LS-jblock,blockSize_);
                    A2AMatrixSet<ComplexF> block(buf_.data(), nExt_ , nStr_ , nt, iblockSize, jblockSize);
                    // std::cout << "  Computing block at (" << iblock << "," << jblock << ")..." << std::endl;

                    LOG(Message) << "Distil matrix block " 
                    << jblock/blockSize_ + nblocki*iblock/blockSize_ + 1 
                    << "/" << nblocki*nblockj << " [" << iblock <<" .. " 
                    << iblock+iblockSize-1 << ", " << jblock <<" .. " << jblock+jblockSize-1 << "]" 
                    << std::endl;

                    int flops       = 0.0;
                    int bytes       = 0.0;
                    int time_kernel = 0.0;
                    double nodes    = env().getGrid()->NodeCount();

                    // loop over cache_ blocks in the current block
                    for(int icache=0 ; icache<iblockSize ; icache+=cacheSize_)   //set according to cache_ size
                    for(int jcache=0 ; jcache<jblockSize ; jcache+=cacheSize_)
                    {
                        int icacheSize = MIN(iblockSize-icache,cacheSize_);      // icacheSize is the size of the current cache_ block (indexed by ii); N_ii-ii is the size of the remainder cache_ block
                        int jcacheSize = MIN(jblockSize-jcache,cacheSize_);
                        A2AMatrixSet<Complex> blockCache(cache_.data(), nExt_, nStr_, nt, icacheSize, jcacheSize);
                        // std::cout << "      Chunk (" << icache << "," << jcache << ")..." << std::endl;

                        double timer=0;
                        startTimer("kernel");
                        A2Autils<FIMPL>::MesonField(blockCache, &left[dtL*dilutionSize_LS+iblock+icache], &right[dtR*dilutionSize_LS+jblock+jcache], gamma_, phase, Tdir, &timer);
                        stopTimer("kernel");

                        time_kernel += timer;

                        flops += vol*(2*8.0+6.0+8.0*nExt_)*icacheSize*jcacheSize*nStr_;
                        bytes += vol*(12.0*sizeof(ComplexD))*icacheSize*jcacheSize
                              +  vol*(2.0*sizeof(ComplexD)*nExt_)*icacheSize*jcacheSize*nStr_;

                        // loop through the cacheblock (inside them) and point blockCache to block
                        startTimer("cache copy");
                        if(par().MesonFieldCase=="phi-phi")
                        {
                            thread_for_collapse( 5, iExt ,nExt_,{
                            for(int iStr=0 ;iStr<nStr_ ; iStr++)
                            for(int t=0 ; t<nt ; t++)
                            for(int iicache=0 ; iicache<icacheSize ; iicache++)
                            for(int jjcache=0;  jjcache<jcacheSize ; jjcache++)
                                block(iExt,iStr,t,icache+iicache,jcache+jjcache) = blockCache(iExt,iStr,t,iicache,jjcache);
                            });
                        }
                        else
                        {
                            thread_for_collapse( 5, iExt ,nExt_,{
                            for(int iStr=0 ;iStr<nStr_ ; iStr++)
                            for(int it=0 ; it<nt_nonzero_ ; it++)  //only wish to copy non-zero timeslices to block (referencing LEFT-rho dilution scheme)
                            for(int iicache=0 ; iicache<icacheSize ; iicache++)
                            for(int jjcache=0;  jjcache<jcacheSize ; jjcache++)
                                block(iExt,iStr,it,icache+iicache,jcache+jjcache) = blockCache(iExt,iStr,lrSourceTimes_.at("left")[dtL][it],iicache,jjcache);
                            });
                        }
                        stopTimer("cache copy");

                        
                    }

                    LOG(Message) << "Kernel perf (flops) " << flops/time_kernel/1.0e3/nodes 
                                << " Gflop/s/node " << std::endl;
                    LOG(Message) << "Kernel perf (read) " << bytes/time_kernel*1.0e6/1024/1024/1024/nodes 
                                << " GB/s/node "  << std::endl;
                    
                // std::cout << " " << bytes << " " << vol  << " " << nExt_ << " " << nStr_ << " " << sizeof(ComplexD) << std::endl;

                    // saving current block to disk
                    startTimer("io");
#ifdef HADRONS_A2AM_PARALLEL_IO
                    //parallel io
                    int inode = env().getGrid()->ThisRank();
                    int nnode = env().getGrid()->RankCount(); 
                    LOG(Message) << "Starting parallel IO. Rank count=" << nnode  << std::endl;
                    env().getGrid()->Barrier();
                    for(int ies=inode ; ies<nExt_*nStr_ ; ies+=nnode){
                        int iExt = ies/nStr_;
                        int iStr = ies%nStr_;
                        if(iblock==0 && jblock==0){              // creates dataset only if it's the first block of the dataset
                            matrixIoTable[iStr + nStr_*iExt].saveBlock(block, iExt , iStr , iblock, jblock, datasetName, cacheSize_);   //set surface chunk size as cacheSize_ (the chunk itself is 3D)
                        }
                        else{
                            matrixIoTable[iStr + nStr_*iExt].saveBlock(block, iExt , iStr , iblock, jblock, datasetName);
                        }
                    }
                    env().getGrid()->Barrier();
#else
                    // serial io
                    LOG(Message) << "Starting serial IO" << std::endl;
                    for(int iExt=0; iExt<nExt_; iExt++)
                    for(int iStr=0; iStr<nStr_; iStr++)
                    {
                        if(iblock==0 && jblock==0){              // creates dataset only if it's the first block of the dataset
                            matrixIoTable[iStr + nStr_*iExt].saveBlock(block, iExt, iStr, iblock, jblock, datasetName, cacheSize_);   //set surface chunk size as cacheSize_ (the chunk itself is 3D)
                        }
                        else{
                            matrixIoTable[iStr + nStr_*iExt].saveBlock(block, iExt, iStr, iblock, jblock, datasetName);
                        }
                    }
#endif
                    stopTimer("io");
                }
            }
        }
    }
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DistilMesonField_hpp_
