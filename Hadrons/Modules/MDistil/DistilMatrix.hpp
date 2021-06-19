#ifndef Distil_matrix_hpp_
#define Distil_matrix_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>

#ifndef HADRONS_DISTIL_IO_TYPE
#define HADRONS_DISTIL_IO_TYPE ComplexF
#endif

#define DISTIL_MATRIX_NAME      "DistilMesonField"
#define METADATA_NAME           "Metadata"
#define DILUTION_METADATA_NAME  "DilutionSchemes"

#define START_TIMER(name) if (tarray) tarray->startTimer(name)
#define STOP_TIMER(name)  if (tarray) tarray->stopTimer(name)
#define GET_TIMER(name)   ((tarray != nullptr) ? tarray->getDTimer(name) : 0.)

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

using DilutionMap  = std::array<std::vector<std::vector<unsigned int>>,3>;
enum Side {left = 0, right = 1};
const std::vector<Side> sides =  {Side::left,Side::right};  //to facilitate iteration

// metadata serialiser class
template <typename FImpl>
class DistilMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldMetadata,
                                    unsigned int,               Nt,
                                    unsigned int,               Nvec,
                                    std::vector<RealF>,         Momentum,
                                    Gamma::Algebra,             Operator,               // just gamma matrices for now, but could turn into more general operators in the future
                                    std::vector<unsigned int>,  NoisePair,
                                    std::string,                MesonFieldType,
                                    std::vector<std::string>,   NoiseHashLeft,
                                    std::vector<std::string>,   NoiseHashRight)
};

//metadata io class to deal with 2d ragged arrays (possibly remove after Mike's serialisation rewrite)
class DistilMetadataIo
{
private:
    std::string fileName_, metadataName_;
public:
    DistilMetadataIo(std::string filename, std::string metadataname):fileName_(filename) , metadataName_(metadataname) {}
    template <typename T>
    void write2dMetadata(const std::string name, const std::vector<std::vector<T>>& data, const std::string newgroupname="")
    {
#ifdef HAVE_HDF5
        // variable-length struct
        typedef struct  {
            size_t len; /* Length of VL data (in base type units) */      
            void *p;    /* Pointer to VL data */        
        } VlStorage;
        
        Hdf5Reader  reader(fileName_, false);
        push(reader, metadataName_);
        H5NS::Group &subgroup = reader.getGroup();

        if(!newgroupname.empty())
        {
            H5NS::Exception::dontPrint();
            try{
                subgroup = subgroup.openGroup(newgroupname);
            } catch (...) {
                subgroup = subgroup.createGroup(newgroupname);
            }
        }

        H5NS::VarLenType vl_type(Hdf5Type<T>::type());
        std::vector<VlStorage> vl_data(data.size());
        for(unsigned int i=0 ; i<data.size() ; i++)
        {
            vl_data.at(i).len = data.at(i).size();
            vl_data.at(i).p = (void*) &data.at(i).front();
        }

        hsize_t         attrDim = data.size();
        H5NS::DataSpace attrSpace(1, &attrDim);
        H5NS::Attribute attr = subgroup.createAttribute(name, vl_type, attrSpace);
        attr.write(vl_type, &vl_data.front());
#else
        HADRONS_ERROR(Implementation, "distillation I/O needs HDF5 library");
#endif
    }
};

//computation class declaration
template <typename FImpl, typename T, typename Tio>
class DmfComputation
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef DistillationNoise<FImpl> DistillationNoise;
    typedef Vector<FermionField> DistilVector;
    typedef typename DistillationNoise::Index Index;
    typedef typename DistillationNoise::LapPack LapPack;
    typedef std::function<std::string(const unsigned int, const unsigned int, const int, const int)>  FilenameFn;
    typedef std::function<DistilMesonFieldMetadata<FImpl>(const unsigned int, const unsigned int, const int, const int)>  MetadataFn;
public:
    long    global_counter = 0;
    double  global_flops = 0.0, global_bytes = 0.0, global_iospeed = 0.0;
private:
    std::map<Side,std::string>          dmfType_;
    GridCartesian*                      g_;
    GridCartesian*                      g3d_;
    ColourVectorField                   evec3d_;
    FermionField                        tmp3d_;
    Vector<Tio>                         bBuf_;
    Vector<T>                           cBuf_;
    const unsigned int                  bSize_;
    const unsigned int                  cSize_;
    const unsigned int                  nt_;
    const unsigned int                  nd_;
    const unsigned int                  next_;
    const unsigned int                  nstr_;
    const bool                          is_exact_;
    std::map<Side, unsigned int>        dil_size_ls_;
    std::map<Side, DistillationNoise&>  noises_;
public:
    DmfComputation(std::map<Side,std::string>   mftype,
                   GridCartesian*               g,
                   GridCartesian*               g3d,
                   DistillationNoise&           nl,
                   DistillationNoise&           nr,
                   const unsigned int           blockSize,
                   const unsigned int           cacheSize,
                   const unsigned int           nt,
                   const unsigned int           next,
                   const unsigned int           nstr,
                   const bool                   is_exact);
    bool isPhi(Side s);
    bool isRho(Side s);
    DilutionMap getMap(Side s);
public:
    void makePhiComponent(FermionField&         phiComponent,
                          DistillationNoise&    n,
                          const unsigned int    n_idx,
                          const unsigned int    iD,
                          PerambTensor&         peramb,
                          LapPack&              epack);
    void makeRhoComponent(FermionField&         rhoComponent,
                          DistillationNoise&    n,
                          const unsigned int    n_idx,
                          const unsigned int    iD);
    void makeDistVecs(std::map<Side, DistilVector&>                 dv,
                      std::vector<unsigned int>                     n_idx,
                      LapPack&                                      epack,
                      std::map<Side, std::vector<unsigned int>>     timeDilSource,
                      std::map<Side, PerambTensor&>                 peramb={});
    void execute(const FilenameFn                               &filenameDmfFn,
                 const MetadataFn                               &metadataDmfFn,
                 std::vector<Gamma::Algebra>                    gamma_,
                 std::map<Side, DistilVector&>                  dv,
                 std::vector<unsigned int>                      n_idx,
                 std::vector<ComplexField>                      ph,
                 std::map<Side, std::vector<unsigned int>>      timeDilSource,
                 TimerArray*                                    tarray);
};

//####################################
//# computation class implementation #
//####################################

template <typename FImpl, typename T, typename Tio>
DmfComputation<FImpl,T,Tio>
::DmfComputation(std::map<Side,std::string>  mftype,
                 GridCartesian*         g,
                 GridCartesian*         g3d,
                 DistillationNoise&     nl,
                 DistillationNoise&     nr,
                 const unsigned int     blockSize,
                 const unsigned int     cacheSize,
                 const unsigned int     nt,
                 const unsigned int     next,
                 const unsigned int     nstr,
                 const bool             is_exact)
: dmfType_(mftype), g_(g), g3d_(g3d), evec3d_(g3d), tmp3d_(g3d)
, nt_(nt) , nd_(g->Nd()), bSize_(blockSize) , cSize_(cacheSize)
, next_(next) , nstr_(nstr) , is_exact_(is_exact)
{
    cBuf_.resize(next_*nstr_*nt_*cSize_*cSize_);
    bBuf_.resize(next_*nstr_*nt_*bSize_*bSize_); //maximum size

    noises_ = std::map<Side, DistillationNoise&> ({{Side::left,nl},{Side::right,nr}});
    dil_size_ls_ = { {Side::left,nl.dilutionSize(Index::l)*nl.dilutionSize(Index::s)} ,
                             {Side::right,nr.dilutionSize(Index::l)*nr.dilutionSize(Index::s)} };
}

template <typename FImpl, typename T, typename Tio>
DilutionMap DmfComputation<FImpl,T,Tio>::getMap(Side s)
{
    DilutionMap m;
    for(auto dil_idx : { Index::t, Index::l, Index::s })
    for(unsigned int it=0 ; it<noises_.at(s).dilutionSize(dil_idx) ; it++)
    {
        std::vector<unsigned int> temp = noises_.at(s).getDilutionPartition(dil_idx,it);
        m[dil_idx].push_back(temp);
    }
    return m;
}

template <typename FImpl, typename T, typename Tio>
bool DmfComputation<FImpl,T,Tio>::isPhi(Side s)
{
    return (dmfType_.at(s)=="phi" ? true : false);
}

template <typename FImpl, typename T, typename Tio>
bool DmfComputation<FImpl,T,Tio>::isRho(Side s)
{
    return (dmfType_.at(s)=="rho" ? true : false);
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makePhiComponent(FermionField&            phiComponent,
                   DistillationNoise&       n,
                   const unsigned int       n_idx,
                   const unsigned int       iD,
                   PerambTensor&            peramb,
                   LapPack&                 epack)
{
    std::array<unsigned int,3> mftype = n.dilutionCoordinates(iD);
    unsigned int dt = mftype[Index::t] , dl = mftype[Index::l] , ds = mftype[Index::s];
    std::vector<int> p_ts = peramb.MetaData.timeSources;
    std::vector<int>::iterator itr_dt = std::find(p_ts.begin(), p_ts.end(), dt);
    unsigned int idt = std::distance(p_ts.begin(), itr_dt); //gets correspondent index of dt in the tensor obj
    const unsigned int nVec = epack.evec.size();
    const unsigned int Ntfirst = g_->LocalStarts()[nd_ - 1];
    const unsigned int Ntlocal = g_->LocalDimensions()[nd_ - 1];
    for (unsigned int t = Ntfirst; t < Ntfirst + Ntlocal; t++)
    {
        tmp3d_ = Zero();
        for (unsigned int k = 0; k < nVec; k++)
        {
            ExtractSliceLocal(evec3d_,epack.evec[k],0,t-Ntfirst,nd_ - 1);
            tmp3d_ += evec3d_ * peramb.tensor(t, k, dl, n_idx, idt, ds);
        }
        InsertSliceLocal(tmp3d_,phiComponent,0,t-Ntfirst,nd_ - 1);
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeRhoComponent(FermionField&        rhoComponent,
                   DistillationNoise&   n,
                   const unsigned int            n_idx,
                   const unsigned int   iD)   
{
    rhoComponent = n.makeSource(iD, n_idx);
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeDistVecs(std::map<Side, DistilVector&>            dv,
          std::vector<unsigned int>                     n_idx,
          LapPack&                                      epack,
          std::map<Side, std::vector<unsigned int>>     timeDilSource,
          std::map<Side, PerambTensor&>                 peramb)
{
    std::map<Side,unsigned int> noise_pair = {{Side::left,n_idx[0]},{Side::right,n_idx[1]}};

    for(auto s : sides)    // computation of phi or rhos
    {
        std::vector<unsigned int> time_sources = timeDilSource.at(s);
        for(unsigned int idt=0 ; idt<time_sources.size() ; idt++)   // compute just time sources specified 
        {
            unsigned int dt = time_sources[idt];
            unsigned int iD_offset = noises_.at(s).dilutionIndex(dt,0,0);    //assuming t is the slowest index
            unsigned int iiD = idt*dil_size_ls_.at(s);  // reduced dilution index
            for(unsigned int iD=iD_offset ; iD<dil_size_ls_.at(s) ; iD++)
            {
                if(isPhi(s))
                {
                    makePhiComponent(dv.at(s)[iiD] , noises_.at(s) , noise_pair.at(s) , iD , peramb.at(s), epack);
                }
                else if(isRho(s))
                {
                    makeRhoComponent(dv.at(s)[iiD] , noises_.at(s) , noise_pair.at(s), iD);
                }
            }
        }
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::execute(const FilenameFn                              &filenameDmfFn,
          const MetadataFn                              &metadataDmfFn,
          std::vector<Gamma::Algebra>                   gamma_,
          std::map<Side, DistilVector&>                 dv,
          std::vector<unsigned int>                     n_idx,
          std::vector<ComplexField>                     ph,
          std::map<Side, std::vector<unsigned int>>     timeDilSource,
          TimerArray*                                   tarray)
{
    std::vector<std::vector<unsigned int>> timeDilutionPairList;
    bool fileIsInit = false;
    const unsigned int vol = g_->_gsites;
    // computing time-dillution blocks and saving to disk
    for (unsigned int idtL=0 ; idtL<timeDilSource.at(Side::left).size() ; idtL++)
    for (unsigned int idtR=0 ; idtR<timeDilSource.at(Side::right).size() ; idtR++)  
    {
        unsigned int dtL = timeDilSource.at(Side::left)[idtL];
        unsigned int dtR = timeDilSource.at(Side::right)[idtR];
        // fetch necessary time slices for this block
        std::map<Side,std::vector<unsigned int>> p = { {Side::left,{}} , {Side::right,{}}};
        for(auto s : sides)
        {
            if(isPhi(s))
            {
                p.at(s).resize(nt_);
                //phi: filling with all time slices -> intersects with non-empty
                std::iota(std::begin(p.at(s)), std::end(p.at(s)), 0); 
            }
            else
            {
                p.at(s) =  noises_.at(s).timeSlices(s==Side::left ? dtL : dtR);
            }
        }
        std::vector<unsigned int> stInter;
        std::set_intersection(p.at(Side::left).begin(), p.at(Side::left).end(), 
                            p.at(Side::right).begin(), p.at(Side::right).end(),
                            std::back_inserter(stInter));

        const int nt_sparse = stInter.size();
        bBuf_.resize(next_*nstr_*nt_sparse*bSize_*bSize_);

        if( !stInter.empty() ) // only execute case when partitions have at least one time slice in common
        {
            timeDilutionPairList.push_back({dtL,dtR});
            LOG(Message) << "------------------------ " << dtL << "-" << dtR << " ------------------------" << std::endl; 
            LOG(Message) << "Saving time slices : " << MDistil::timeslicesDump(stInter) << std::endl;
            LOG(Message) << "Time extension in file : " << nt_sparse << std::endl;
            std::string datasetName = std::to_string(dtL)+"-"+std::to_string(dtR);

            unsigned int nblocki = dil_size_ls_.at(Side::left)/bSize_ + (((dil_size_ls_.at(Side::left) % bSize_) != 0) ? 1 : 0);
            unsigned int nblockj = dil_size_ls_.at(Side::right)/bSize_ + (((dil_size_ls_.at(Side::right) % bSize_) != 0) ? 1 : 0);

            // loop over blocks in the current time-dilution block
            for(unsigned int iblock=0 ; iblock<dil_size_ls_.at(Side::left) ; iblock+=bSize_) //set according to memory size
            for(unsigned int jblock=0 ; jblock<dil_size_ls_.at(Side::right) ; jblock+=bSize_)
            {
                unsigned int iblockSize = MIN(dil_size_ls_.at(Side::left)-iblock,bSize_);    // iblockSize is the size of the current block (indexed by i); N_i-i is the size of the eventual remainder block
                unsigned int jblockSize = MIN(dil_size_ls_.at(Side::right)-jblock,bSize_);
                A2AMatrixSet<Tio> block(bBuf_.data(), next_ , nstr_ , nt_sparse, iblockSize, jblockSize);

                LOG(Message) << "Distil matrix block " 
                << jblock/bSize_ + nblocki*iblock/bSize_ + 1 
                << "/" << nblocki*nblockj << " [" << iblock << " .. " 
                << iblock+iblockSize-1 << ", " << jblock << " .. " << jblock+jblockSize-1 << "]" 
                << std::endl;

                double flops        = 0.0;
                double bytes        = 0.0;
                double time_kernel  = 0.0;
                double nodes        = g_->NodeCount();

                // loop over cache_ blocks in the current block
                for(unsigned int icache=0 ; icache<iblockSize ; icache+=cSize_)   //set according to cache_ size
                for(unsigned int jcache=0 ; jcache<jblockSize ; jcache+=cSize_)
                {
                    unsigned int icacheSize = MIN(iblockSize-icache,cSize_);      // icacheSize is the size of the current cache_ block (indexed by ii); N_ii-ii is the size of the remainder cache_ block
                    unsigned int jcacheSize = MIN(jblockSize-jcache,cSize_);
                    A2AMatrixSet<T> blockCache(cBuf_.data(), next_, nstr_, nt_, icacheSize, jcacheSize);
                    
                    // temporary log messages
                    // LOG(Message) << "Block size (one be saved to disk, sparse time dimension): "    << nt_sparse*iblockSize*jblockSize*sizeof(Tio)/1024. << "KB/momentum/gamma" << std::endl;
                    // LOG(Message) << "Cache block size (non-sparse time dimension): "                << nt_*icacheSize*jcacheSize*sizeof(T)/1024. << "KB/momentum/gamma" << std::endl;

                    double timer = 0.0;
                    START_TIMER("kernel");
                    // as dictated by DistillationNoise, dt must be the slowest index for this to work; otherwise will have to compute l/r block at each contraction)
                    unsigned int iiDL = idtL*dil_size_ls_.at(Side::left) , iiDR = idtR*dil_size_ls_.at(Side::right);
                    A2Autils<FImpl>::MesonField(blockCache, &dv.at(Side::left)[iiDL+iblock+icache], &dv.at(Side::right)[iiDR+jblock+jcache], gamma_, ph, nd_ - 1, &timer);
                    STOP_TIMER("kernel");
                    // std::cout << "intermediate kernel time: " << timer << " ; " << "grid kernel time: " << GET_TIMER("kernel") << std::endl;
                    // std::cin.get();
                    time_kernel += timer;

                    // nExt is currently # of momenta , nStr is # of gamma matrices
                    flops += vol*(2*8.0+6.0+8.0*next_)*icacheSize*jcacheSize*nstr_;
                    bytes += vol*(12.0*sizeof(T))*icacheSize*jcacheSize
                            +  vol*(2.0*sizeof(T)*next_)*icacheSize*jcacheSize*nstr_;

                    // loop through the cacheblock (inside them) and point blockCache to block
                    START_TIMER("cache copy");
                    unsigned int stSize = stInter.size();
                    thread_for_collapse(5,iext,next_,{
                    for(unsigned int istr=0;istr<nstr_;istr++)
                    for(unsigned int it=0;it<stSize;it++)
                    for(unsigned int iicache=0;iicache<icacheSize;iicache++)
                    for(unsigned int jjcache=0;jjcache<jcacheSize;jjcache++)
                    {
                        block(iext,istr,it,icache+iicache,jcache+jjcache)=blockCache(iext,istr,stInter[it],iicache,jjcache);
                    }
                    });
                    STOP_TIMER("cache copy");
                }

                LOG(Message) << "Kernel perf (flops) " << flops/time_kernel/1.0e3/nodes 
                            << " Gflop/s/node " << std::endl;
                LOG(Message) << "Kernel perf (read) " << bytes/time_kernel*0.000931322574615478515625/nodes //  *1.0e6/1024/1024/1024/nodes
                            << " GB/s/node "  << std::endl;
                global_counter++;
                global_flops += flops/time_kernel/1.0e3/nodes ;
                global_bytes += bytes/time_kernel*0.000931322574615478515625/nodes ; // 1.0e6/1024/1024/1024/nodes

                // saving current block to disk
                double ioTime = -GET_TIMER("IO: write block");
                START_TIMER("IO: total");
#ifdef HADRONS_A2AM_PARALLEL_IO
                //parallel io
                unsigned int inode = g_->ThisRank();
                unsigned int nnode = g_->RankCount(); 
                LOG(Message) << "Starting parallel IO. Rank count=" << nnode  << std::endl;
                g_->Barrier();
                for(unsigned int ies=inode ; ies<next_*nstr_ ; ies+=nnode){
                    unsigned int iext = ies/nstr_;
                    unsigned int istr = ies%nstr_;
                    // metadata;
                    DistilMesonFieldMetadata<FImpl> md = metadataDmfFn(iext,istr,n_idx[0],n_idx[1]);                    
                    A2AMatrixIo<HADRONS_DISTIL_IO_TYPE> matrixIo(filenameDmfFn(iext,istr,n_idx[0],n_idx[1]), DISTIL_MATRIX_NAME, nt_sparse, dil_size_ls_.at(Side::left), dil_size_ls_.at(Side::right));
                    START_TIMER("IO: write block");
                    if(iblock==0 && jblock==0)  
                    {             
                        if( (dtL==timeDilSource.at(Side::left)[0]) && (dtR==timeDilSource.at(Side::right)[0]) )     //execute this once per block
                        {
                            START_TIMER("IO: file creation");
                            matrixIo.initFile(md);
                            STOP_TIMER("IO: file creation");
                        }
                        matrixIo.saveBlock(block, iext , istr , iblock, jblock, datasetName, stInter, cSize_);   //sets 2D chunk size and creates dataset
                    }
                    else{
                        matrixIo.saveBlock(block, iext , istr , iblock, jblock, datasetName);
                    }
                    STOP_TIMER("IO: write block");
                }
                g_->Barrier();
#endif
                STOP_TIMER("IO: total");
                ioTime    += GET_TIMER("IO: write block");
                unsigned int bytesBlockSize  = static_cast<double>(next_*nstr_*nt_sparse*iblockSize*jblockSize*sizeof(Tio));
                double iospeed = bytesBlockSize/ioTime*0.95367431640625;     // 1.0e6/1024/1024
                LOG(Message)    << "HDF5 IO done " << sizeString(bytesBlockSize) << " in "
                                << ioTime  << " us (" << iospeed << " MB/s)" << std::endl;
                global_iospeed += iospeed;
            }
        }
    }
    //saving dilution schemes and time source pairs
    if(g_->IsBoss())
    {
        START_TIMER("IO: total");
        for(unsigned int iext=0 ; iext<next_ ; iext++)
        for(unsigned int istr=0 ; istr<nstr_ ; istr++)
        {
            DistilMetadataIo mdIo(filenameDmfFn(iext,istr,n_idx[0],n_idx[1]),
                 std::string(DISTIL_MATRIX_NAME) + "/" + std::string(METADATA_NAME) );
            mdIo.write2dMetadata("TimeSourcePairs", timeDilutionPairList);
            //  dilution schemes (2d ragged metadata)
            DilutionMap lmap = getMap(Side::left);
            DilutionMap rmap = getMap(Side::right);
            mdIo.write2dMetadata("TimeDilutionLeft" , lmap[Index::t], DILUTION_METADATA_NAME);
            mdIo.write2dMetadata("TimeDilutionRight", rmap[Index::t], DILUTION_METADATA_NAME);
            mdIo.write2dMetadata("LapDilutionLeft"  , lmap[Index::l], DILUTION_METADATA_NAME);
            mdIo.write2dMetadata("LapDilutionRight" , rmap[Index::l], DILUTION_METADATA_NAME);
            mdIo.write2dMetadata("SpinDilutionLeft" , lmap[Index::s], DILUTION_METADATA_NAME);
            mdIo.write2dMetadata("SpinDilutionRight", rmap[Index::s], DILUTION_METADATA_NAME);
        }
        STOP_TIMER("IO: total");
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Distil_matrix_hpp_
