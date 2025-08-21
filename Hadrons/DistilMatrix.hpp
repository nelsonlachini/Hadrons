/*
 * DistilMatrix.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Nelson Lachini <nelson.lachini@ed.ac.uk> 
 * Author: Simon Bürger <simon.buerger@rwth-aachen.de>
 * Author: nelsonlachini <nelsonlachini@gmail.com>
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

#ifndef Distil_matrix_hpp_
#define Distil_matrix_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Modules/MDistil/DistilUtils.hpp>
#include <Hadrons/DistillationVectors.hpp>

#ifndef HADRONS_DISTIL_IO_TYPE
#define HADRONS_DISTIL_IO_TYPE ComplexF
#endif

#ifndef HADRONS_DISTIL_TYPE
#define HADRONS_DISTIL_TYPE ComplexD
#endif

// number of distil vector time-dilution components held in memory at the same time
#ifndef DISTILVECTOR_TIME_BATCH_SIZE
#define DISTILVECTOR_TIME_BATCH_SIZE 1
#endif

#define DISTIL_MATRIX_NAME      "DistilMesonField"
#define METADATA_NAME           "Metadata"
#define DILUTION_METADATA_NAME  "DilutionSchemes"

#define START_TIMER(name) if (tarray) tarray->startTimer(name)
#define STOP_TIMER(name)  if (tarray) tarray->stopTimer(name)
#define GET_TIMER(name)   ((tarray != nullptr) ? tarray->getDTimer(name) : 0.)

#define HADRONS_DISTIL_PARALLEL_IO

typedef unsigned int uint;

BEGIN_HADRONS_NAMESPACE

// general distil matrix set for io buffers
// Dimensions:
//   0 - ext/str - external field (momentum, EM field, ...) and operator (gamma matrix, ...) structure
//   1 - t   - timeslice
//   2 - i   - left  spin-Laplacian mode index
//   3 - j   - right spin-Laplacian mode index
template <typename T>
using DistilMatrixSetIo = Eigen::TensorMap<Eigen::Tensor<T, 4, Eigen::RowMajor>>;

// Distil matrix set on a single timeslice for io buffers 
// Dimensions:
//   0 - ext/str - external field (momentum, EM field, ...) and spin structure
//   1 - i   - left  spin-Laplacian mode index
//   2 - j   - right spin-Laplacian mode index
template <typename T>
using DistilMatrixSetTimeSliceIo = Eigen::TensorMap<Eigen::Tensor<T, 3, Eigen::RowMajor>>;

// Distil matrix set passed to A2AUtils::MesonField() kernel ; cache object
template <typename T>
using DistilMatrixSetCache = A2AMatrixSet<T>;

// Distil meson field (on single timeslice) ; used for loading a single time dilution block 
// Dimensions
//   0 : i   - left  spin-Laplacian mode index
//   1 : j   - right spin-Laplacian mode index
template <typename T>
using DistilMesonFieldMatrix = Eigen::Matrix<T, -1, -1, Eigen::RowMajor>;

using DilutionMap  = std::array<std::vector<std::vector<uint>>,3>;

enum Side {left = 0, right = 1};
const std::vector<Side> sides =  {Side::left,Side::right};
const std::map<Side, std::string> sideToStr = {{Side::left,"left"}, {Side::right,"right"}};

// metadata serialiser class
template <typename FImpl>
class DistilMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilMesonFieldMetadata,
                                    uint,                               Nt,
                                    uint,                               Nvec,
                                    std::vector<RealF>,                         Momentum,
                                    Gamma::Algebra,                             Operator,               // potentially more general operators in the future
                                    std::vector<uint>,                  NoisePair,
                                    std::string,                                MesonFieldType,
                                    std::string,                                NoiseHashLeft,
                                    std::string,                                NoiseHashRight,
                                    std::vector<std::vector<uint>>,     TimeDilutionLeft,
                                    std::vector<std::vector<uint>>,     TimeDilutionRight,
                                    std::vector<std::vector<uint>>,     LapDilutionLeft,
                                    std::vector<std::vector<uint>>,     LapDilutionRight,
                                    std::vector<std::vector<uint>>,     SpinDilutionLeft,
                                    std::vector<std::vector<uint>>,     SpinDilutionRight,
                                    std::string,                                RelativeSide,
                                    )
};

/******************************************************************************
 *                  Class to handle Distil matrix block HDF5 I/O              *
 ******************************************************************************/
template <typename T>
class DistilMatrixIo
{
public:
    // constructors
    DistilMatrixIo(void) = default;
    DistilMatrixIo(std::string filename, std::string dataname, const uint nt=0, const uint ni = 0,
                const uint nj = 0);
    // destructor
    ~DistilMatrixIo(void) = default;
    // access
    size_t       getSize(void) const;

    // block I/O
    template <typename MetadataType>
    void initFile(const MetadataType &d, const bool verbose=false);

    void saveBlock(const T *data, const uint i, const uint j,
                   const uint blockSizei, const uint blockSizej, std::string t_name, std::string datasetName, const uint chunkSize);

    void saveBlock(const DistilMatrixSetTimeSliceIo<T> &m, const uint iextstr,
                               const uint i, const uint j, std::string datasetName,
                               const uint t, const uint chunkSize);

    template <template <class> class Vec, typename VecT>    // compatibility with A2A
    void load(Vec<VecT> &v, const uint t, const std::string dataset_name, double *tRead = nullptr, GridBase *grid = nullptr);
    template <typename Mat>
    void load(Mat &v, const uint t, const std::string dataset_name, double *tRead = nullptr, GridBase *grid = nullptr);
private:
    std::string  filename_{""}, dataname_{""};
    uint nt_{0}, ni_{0}, nj_{0};
};

// implementation /////////////////////////////////////////////////////////////////
template <typename T>
DistilMatrixIo<T>::DistilMatrixIo(std::string filename, std::string dataname, 
                            const uint nt, const uint ni,
                            const uint nj)
: filename_(filename), dataname_(dataname), ni_(ni), nj_(nj), nt_(nt)
{}

template <typename T>
template <typename MetadataType>
void DistilMatrixIo<T>::initFile(const MetadataType &d, const bool verbose)
{
#ifdef HAVE_HDF5
    if(!std::filesystem::exists(filename_))
    {
        if(verbose)
        {
            LOG(Message) << "Creating files " << filename_ << " , ..." << std::endl;
        }
        Hdf5Writer writer(filename_);
        push(writer, dataname_);
        write(writer, "Metadata", d);
    }
    else if(verbose)
    {
        LOG(Message) << "WARNING: appending to files " << filename_ << " , ... generated on past executions (suppressing further warnings)" << std::endl;
    }
#else
    HADRONS_ERROR(Implementation, "distil matrix I/O needs HDF5 library");
#endif
    }

template <typename T>
size_t DistilMatrixIo<T>::getSize(void) const
{
    return ni_*nj_*sizeof(T);
}

template <typename T>
void DistilMatrixIo<T>::saveBlock(const T *data, 
                               const uint i, 
                               const uint j,
                               const uint blockSizei,
                               const uint blockSizej,
                               std::string t_name,
                               std::string datasetName,
                               const uint chunkSize)
{
#ifdef HAVE_HDF5
    H5NS::H5File file(filename_,H5F_ACC_RDWR);
    H5NS::Group rootgroup = file.openGroup(DISTIL_MATRIX_NAME);

    H5NS::Group tgroup;
    if ( H5Lexists( rootgroup.getId(), t_name.c_str(), H5P_DEFAULT ) > 0 )
    {
        tgroup = rootgroup.openGroup(t_name);
    }
    else
    {
        tgroup = rootgroup.createGroup(t_name);
    }

    H5NS::DataSet        dataset;
    if ( H5Lexists( tgroup.getId(), datasetName.c_str(), H5P_DEFAULT ) > 0 )
    {
        dataset = tgroup.openDataSet(datasetName);
    }
    else
    {
        H5NS::DataSpace      dataspace;
        std::vector<hsize_t>    dim = {static_cast<hsize_t>(ni_), 
                                    static_cast<hsize_t>(nj_)},
                                chunk = {static_cast<hsize_t>(chunkSize), 
                                    static_cast<hsize_t>(chunkSize)};

        dataspace.setExtentSimple(dim.size(), dim.data());
        H5NS::DSetCreatPropList     plist;
        plist.setChunk(chunk.size(), chunk.data());
        plist.setFletcher32();
        dataset = tgroup.createDataSet(datasetName, Hdf5Type<T>::type(), dataspace, plist);
    }

    std::vector<hsize_t> count = {blockSizei, blockSizej},
                         offset = {static_cast<hsize_t>(i),
                                   static_cast<hsize_t>(j)},
                         stride = {1, 1},
                         block  = {1, 1}; 
    H5NS::DataSpace      memspace(count.size(), count.data()), dataspace;
    dataspace   = dataset.getSpace();
    // std::cout << std::endl << "H5 dims= " << ni_ << " " << nj_ << std::endl;
    // std::cout << "H5 chunks= " << chunkSize_i << " " << chunkSize_j << std::endl;
    // std::cout << "H5 counts= " << blockSizei << " " << blockSizej << std::endl;
    // std::cout << "H5 offsets= " << i << " " << j << std::endl;

    dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                              stride.data(), block.data());
    dataset.write(data, Hdf5Type<T>::type(), memspace, dataspace);

#else
    HADRONS_ERROR(Implementation, "distil matrix I/O needs HDF5 library");
#endif
}

template <typename T>
void DistilMatrixIo<T>::saveBlock(const DistilMatrixSetTimeSliceIo<T> &m,
                            //    const uint ext, const uint str,
                               const uint iextstr,  //local
                               const uint i, const uint j, std::string datasetName,
                               const uint t, const uint chunkSize)
{
    uint blockSizei = m.dimension(1);
    uint blockSizej = m.dimension(2);
    // uint nextstr    = m.dimension(1);
    size_t       offset     = (iextstr*nt_ + t)*blockSizei*blockSizej;

    std::string t_name = std::to_string(t);

    saveBlock(m.data() + offset, i, j, blockSizei, blockSizej, t_name, datasetName, chunkSize);
}

template <typename T>
template <template <class> class Vec, typename VecT>
void DistilMatrixIo<T>::load(Vec<VecT> &v, const uint t, const std::string dataset_name, double *tRead, GridBase *grid)
{
#ifdef HAVE_HDF5
    std::vector<hsize_t> hdim;
    H5NS::DataSet        dataset;
    H5NS::DataSpace      dataspace;
    H5NS::CompType       datatype;

    if (!(grid) or grid->IsBoss())
    {
        Hdf5Reader reader(filename_);
        push(reader, dataname_);
        auto &root_group = reader.getGroup();
        std::string t_name = std::to_string(t);
        auto t_group = root_group.openGroup(t_name);
        dataset = t_group.openDataSet(dataset_name);
        datatype = dataset.getCompType();
        dataspace = dataset.getSpace();
        hdim.resize(dataspace.getSimpleExtentNdims());
        dataspace.getSimpleExtentDims(hdim.data());
        if ((ni_ * nj_ != 0) and
            ((hdim[0] != ni_) or (hdim[1] != nj_)))
        {
            HADRONS_ERROR(Size, "distil matrix size mismatch (got "
                + std::to_string(hdim[0]) + "x" + std::to_string(hdim[1]) + ", expected "
                + std::to_string(ni_) + "x" + std::to_string(nj_));
        }
        ni_ = hdim[0];
        nj_ = hdim[1];
    }
    if (grid)
    {
        grid->Broadcast(grid->BossRank(), &ni_, sizeof(uint));
        grid->Broadcast(grid->BossRank(), &nj_, sizeof(uint));
    }

    DistilMesonFieldMatrix<T>         buf(ni_, nj_);
    int broadcastSize =  sizeof(T) * buf.size();
    std::vector<hsize_t> count    = {static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_)},
                         stride   = {1, 1},
                         block    = {1, 1},
                         memCount = {static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_)};
    H5NS::DataSpace      memspace(memCount.size(), memCount.data());

    std::cout << "Loading timeslice";
    std::cout.flush();
    *tRead = 0.;

    std::vector<hsize_t> offset = {0, 0};
    
    if (!(grid) or grid->IsBoss())
    {
        dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                    stride.data(), block.data());
    }
    if (tRead) *tRead -= usecond();
    if (!(grid) or grid->IsBoss())
    {
        dataset.read(buf.data(), datatype, memspace, dataspace);
    }
    if (grid)
    {
        grid->Broadcast(grid->BossRank(), buf.data(), broadcastSize);
    }
    if (tRead) *tRead += usecond();
    v[0] = buf.template cast<VecT>();
    std::cout << std::endl;
#else
    HADRONS_ERROR(Implementation, "distil matrix I/O needs HDF5 library");
#endif
}

template <typename T>
template <typename Mat>
void DistilMatrixIo<T>::load(Mat &v, const uint t, const std::string dataset_name, double *tRead, GridBase *grid)
{
#ifdef HAVE_HDF5
    std::vector<hsize_t> hdim;
    H5NS::DataSet        dataset;
    H5NS::DataSpace      dataspace;
    H5NS::CompType       datatype;

    Hdf5Reader reader(filename_);
    push(reader, dataname_);
    auto &root_group = reader.getGroup();
    std::string t_name = std::to_string(t);
    auto t_group = root_group.openGroup(t_name);
    dataset = t_group.openDataSet(dataset_name);
    datatype = dataset.getCompType();
    dataspace = dataset.getSpace();
    hdim.resize(dataspace.getSimpleExtentNdims());
    dataspace.getSimpleExtentDims(hdim.data());
    if ((ni_ * nj_ != 0) and
        ((hdim[0] != ni_) or (hdim[1] != nj_)))
    {
        HADRONS_ERROR(Size, "distil matrix size mismatch (got "
            + std::to_string(hdim[0]) + "x" + std::to_string(hdim[1]) + ", expected "
            + std::to_string(ni_) + "x" + std::to_string(nj_));
    }
    ni_ = hdim[0];
    nj_ = hdim[1];

    DistilMesonFieldMatrix<T>         buf(ni_, nj_);
    int broadcastSize =  sizeof(T) * buf.size();
    std::vector<hsize_t> count    = {static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_)},
                         stride   = {1, 1},
                         block    = {1, 1},
                         memCount = {static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_)};
    H5NS::DataSpace      memspace(memCount.size(), memCount.data());

    // LOG(Message) << "Loading timeslice " << std::to_string(t) << " block " << dataset_name << std::endl;
    std::cout.flush();
    *tRead = 0.;

    std::vector<hsize_t> offset = {0, 0};

    dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                stride.data(), block.data());

    if (tRead) *tRead -= usecond();

    dataset.read(buf.data(), datatype, memspace, dataspace);
    if (tRead) *tRead += usecond();
    // v = static_cast<Mat>(buf);
    // v = buf.template cast<VecT>();
    v = buf ;       //bring this up in meeting
#else
    HADRONS_ERROR(Implementation, "distil matrix I/O needs HDF5 library");
#endif
}

//####################################
//# computation class declaration    #
//####################################

template <typename FImpl, typename T, typename Tio>
class DmfComputation
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef Hadrons::DistillationNoise<FImpl> DistillationNoise;
    typedef std::vector<FermionField> DistilVector;
    typedef typename DistillationNoise::Index Index;
    typedef typename DistillationNoise::LapPack LapPack;
    typedef std::function<std::string(const uint, const uint, const int, const int)>  TrajNumberFn;
    typedef std::function<std::string(const uint, const uint, const int, const int)>  FilenameFn;
    typedef std::function<DistilMesonFieldMetadata<FImpl>(const uint, const uint, const int, const int, DilutionMap, DilutionMap)>  MetadataFn;
public:
    long    blockCounter_ = 0;
    double  blockFlops_ = .0, blockBytes_ = .0, blockIoSpeed_ = .0;
private:
    std::map<Side,std::string>          dmfType_;
    GridCartesian*                      g_;
    GridCartesian*                      g3d_;
    ColourVectorField                   evec3d_;
    FermionField                        tmp3d_, tmp4d_;
    const uint                          blockSize_; //eventually turns into io chunk size
    const uint                          cacheSize_;
    const uint                          traj_;
    const uint                          nt_, nd_, nExt_, nStr_;
    const bool                          isExact_;
    // bool                                isInitFile_=false;
    std::map<Side, uint>                dilSizeLS_;
    std::map<Side, DistillationNoise&>  distilNoise_;
    const uint                          dvBatchSize_ = DISTILVECTOR_TIME_BATCH_SIZE;
    std::map<Side, std::string>         vectorStem_={{Side::left,""},{Side::right,""}};
public:
    DmfComputation(std::map<Side,std::string>   mf_type,
                   GridCartesian*               g,
                   GridCartesian*               g3d,
                   DistillationNoise&           nl,
                   DistillationNoise&           nr,
                   const uint                   block_size,
                   const uint                   cache_size,
                   const uint                   nt,
                   const uint                   n_ext,
                   const uint                   n_str,
                   const bool                   is_exact,
                   const uint                   traj,
                   const std::string            left_vector_stem="",
                   const std::string            right_vector_stem="");
    bool isPhi(Side s);
    bool isRho(Side s);
    DilutionMap fetchDilutionMap(Side s);
private:
    void makePhiComponent(FermionField&             phi_component,
                          DistillationNoise&        n,
                          const uint                n_idx,
                          const uint                D,
                          MDistil::PerambTensor&    peramb,
                          LapPack&                  epack,
                          TimerArray*                                 tarray);
    void loadPhiComponent(FermionField&             phi_component,
                          DistillationNoise&        n,
                          const uint                n_idx,
                          const uint                D,
                          std::string               vector_stem,
                          LapPack&                  epack,
                          TimerArray*                                 tarray);

    void makeRhoComponent(FermionField&         rho_component,
                          DistillationNoise&    n,
                          const uint            n_idx,
                          const uint            D,
                          TimerArray*                                 tarray);
    void makeDvLapSpinComponent(FermionField&                           component,
                                std::map<Side, uint>                    n_idx,
                                LapPack&                                epack,
                                Side                                    s,
                                uint                                    D,
                                std::map<Side, MDistil::PerambTensor&>  peramb,
                                TimerArray*                                 tarray);
    void makeDvLapSpinBlock(std::map<Side, DistilVector&>                       dv,
                                      std::map<Side, uint>                      n_idx,
                                      LapPack&                                  epack,
                                      Side                                      s,
                                      uint                                      dt,
                                      uint                                      iibatch,
                                      TimerArray*                                 tarray,
                                      std::map<Side, MDistil::PerambTensor&>    peramb={});
    void makeDvLapSpinCacheBlock(std::map<Side, DistilVector&>          dv,
                               uint                                     dv_idx_offset,
                               std::map<Side, uint>                     n_idx,
                               LapPack&                                 epack,
                               Side                                     s,
                               std::map<Side, MDistil::PerambTensor&>   peramb,
                               TimerArray*                                 tarray);
    void makeDvLapSpinBatch(std::map<Side, DistilVector&>               dv,
                            std::map<Side, uint>                    n_idx,
                            LapPack&                                epack,
                            Side                                    s,
                            std::vector<uint>                       dt_list,
                            std::map<Side, MDistil::PerambTensor&>  peramb,
                            TimerArray*                                 tarray);
    std::vector<uint> fetchDvBatchIdxs(uint                                ibatch,
                                       std::map<Side, std::vector<uint>>   time_dil_sources,
                                       Side                                s,
                                       uint                                shift=0);
    void makeRelativePhiComponent(FermionField&                         phi_component,
                                DistillationNoise&                      n,
                                const uint                              n_idx,
                                const uint                              D,
                                const uint                              delta_t,
                                std::map<Side, MDistil::PerambTensor&>  peramb,
                                Side                                    s,
                                LapPack&                                epack,
                                TimerArray*                                 tarray);
    void makeRelativeRhoComponent(FermionField&         rho_component,
                            DistillationNoise&          n,
                            const uint                  n_idx,
                            const uint                  D,
                            LapPack&                    epack,
                            TimerArray*                                 tarray);
    void makeRelativeDvLapSpinBlock(std::map<Side, DistilVector&>       dv,
                               std::vector<uint>                        dt_list,
                               std::map<Side, uint>                     n_idx,
                               LapPack&                                 epack,
                               Side                                     s,
                               const uint                               delta_t,
                               std::map<Side, MDistil::PerambTensor&>   peramb,
                               TimerArray*                                 tarray);
public:
    void executeRelative(const FilenameFn&                              filenameDmfFn,
                       const MetadataFn&                                metadataDmfFn,
                       Vector<Tio>&                                     bBuf,
                       Vector<T>&                                       cBuf,
                       std::vector<Gamma::Algebra>                      gamma,
                       std::map<Side, DistilVector&>                    dv,
                       std::map<Side, uint>                             n_idx,
                       std::vector<ComplexField>                        ph,
                       std::map<Side, std::vector<uint>>                time_dil_source,
                       LapPack&                                         epack,
                       TimerArray*                                      tarray,
                       Side                                             relative_side,
                       std::vector<uint>                                delta_t_list,
                       std::map<Side, MDistil::PerambTensor&>           peramb={});
    void executeFixed(const FilenameFn&                         filenameDmfFn,
                 const MetadataFn&                              metadataDmfFn,
                 Vector<Tio>&                                   bBuf,
                 Vector<T>&                                     cBuf,
                 std::vector<Gamma::Algebra>                    gamma,
                 std::map<Side, DistilVector&>                  dv,
                 std::map<Side, uint>                           n_idx,
                 std::vector<ComplexField>                      ph,
                 std::map<Side, std::vector<uint>>              time_dil_source,
                 LapPack&                                       epack,
                 TimerArray*                                    tarray,
                 bool                                           only_diag,
                 const uint                                     diag_shift,
                 std::map<Side, MDistil::PerambTensor&>         peramb={});
};

//####################################
//# computation class implementation #
//####################################

template <typename FImpl, typename T, typename Tio>
DmfComputation<FImpl,T,Tio>
::DmfComputation(std::map<Side,std::string>     mf_type,
                 GridCartesian*                 g,
                 GridCartesian*                 g3d,
                 DistillationNoise&             nl,
                 DistillationNoise&             nr,
                 const uint                     block_size,
                 const uint                     cache_size,
                 const uint                     nt,
                 const uint                     n_ext,
                 const uint                     n_str,
                 const bool                     is_exact,
                 const uint                     traj,
                 const std::string              left_vector_stem,
                 const std::string              right_vector_stem)
: dmfType_(mf_type), g_(g), g3d_(g3d), evec3d_(g3d), tmp3d_(g3d), tmp4d_(g)
, nt_(nt) , nd_(g->Nd()), blockSize_(block_size) , cacheSize_(cache_size)
, nExt_(n_ext) , nStr_(n_str) , isExact_(is_exact), traj_(traj)
{
    distilNoise_ = std::map<Side, DistillationNoise&> ({{Side::left,nl},{Side::right,nr}});
    dilSizeLS_ = { {Side::left,nl.dilutionSize(Index::l)*nl.dilutionSize(Index::s)} ,
                             {Side::right,nr.dilutionSize(Index::l)*nr.dilutionSize(Index::s)} };

    vectorStem_ = { {Side::left,left_vector_stem} , {Side::right,right_vector_stem}};
}

template <typename FImpl, typename T, typename Tio>
DilutionMap DmfComputation<FImpl,T,Tio>::fetchDilutionMap(Side s)
{
    DilutionMap m;
    for(Index dil_idx : { Index::t, Index::l, Index::s })
    for(uint it=0 ; it<distilNoise_.at(s).dilutionSize(dil_idx) ; it++)
    {
        std::vector<uint> temp = distilNoise_.at(s).dilutionPartition(dil_idx,it);
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

// fetch time dilution indices (sources) in dv batch ibatch
template <typename FImpl, typename T, typename Tio>
std::vector<uint> DmfComputation<FImpl,T,Tio>
::fetchDvBatchIdxs(uint ibatch, std::map<Side, std::vector<uint>> time_dil_sources, Side s, const uint shift)
{
    std::vector<uint> batch_dt;
    for(uint dt=ibatch*dvBatchSize_; dt<(ibatch+1)*dvBatchSize_; dt++){
        batch_dt.push_back( (time_dil_sources.at(s)[dt]+shift)%distilNoise_.at(s).dilutionSize(Index::t));
    }
    return batch_dt;
}

// DV FIXED (ANCHORED) METHODS
template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makePhiComponent(FermionField&            phi_component,
                   DistillationNoise&       n,
                   const uint               n_idx,
                   const uint               D,
                   MDistil::PerambTensor&   peramb,
                   LapPack&                 epack,
                   TimerArray*                                 tarray)
{
    // LOG(Message) << std::endl << std::endl << "--makePhiComponent(D="<<D<<")" << std::endl;
    std::array<uint,3> d_coor = n.dilutionCoordinates(D);
    uint dt = d_coor[Index::t] , dk = d_coor[Index::l] , ds = d_coor[Index::s];
    std::vector<int> peramb_ts = peramb.MetaData.timeSources;
    std::vector<int>::iterator itr_dt = std::find(peramb_ts.begin(), peramb_ts.end(), dt);
    uint idt = std::distance(peramb_ts.begin(), itr_dt); //gets correspondent index of dt in the tensor obj 
    const uint nVec = epack.evec.size();
    const uint Nt_first = g_->LocalStarts()[nd_ - 1];
    const uint Nt_local = g_->LocalDimensions()[nd_ - 1];
    for (uint t = Nt_first; t < Nt_first + Nt_local; t++)
    {
        tmp3d_ = Zero();
        for (uint k = 0; k < nVec; k++)
        {
            // double start = usecond();
            START_TIMER("makePhiComponent:ExtractSliceLocal");
            ExtractSliceLocal(evec3d_,epack.evec[k],0,t-Nt_first,nd_ - 1);
            STOP_TIMER("makePhiComponent:ExtractSliceLocal");
            // double stop = usecond();
            // LOG(Message) << "ExtractSliceLocal(evec["<<k<<"], slice_lo=0, slice_hi="<<t-Nt_first<<") : " << stop-start << " us" << std::endl;

            // start = usecond();
            START_TIMER("makePhiComponent:tmp3d_+=evec3d_");
            tmp3d_ += evec3d_ * peramb.tensor(t, k, dk, n_idx, idt, ds);
            STOP_TIMER("makePhiComponent:tmp3d_+=evec3d_");
            // stop = usecond();
            // LOG(Message) << "tmp3d_ += evec3d_  (t="<<t<<", dk="<<dk<<", idt="<<idt<<", ds="<<ds<<") : " << stop-start << " us" << std::endl << std::endl;
        }
        // double start = usecond();
        START_TIMER("makePhiComponent:InsertSliceLocal");
        InsertSliceLocal(tmp3d_,phi_component,0,t-Nt_first,nd_ - 1);
        STOP_TIMER("makePhiComponent:InsertSliceLocal");
        // double stop = usecond();
        // LOG(Message) << "InsertSliceLocal(slice_lo=0, slice_hi="<<t-Nt_first<<") : " << stop-start << " us" << std::endl << std::endl;
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::loadPhiComponent(FermionField&            phi_component,
                   DistillationNoise&       n,
                   const uint               n_idx,
                   const uint               D,
                   std::string              vector_stem,
                   LapPack&                 epack,
                   TimerArray*                                 tarray)
{
    const uint nVec = epack.evec.size();
    START_TIMER("IO: total");
    START_TIMER("IO: read");
    DistillationVectorsIo::readComponent(phi_component, vector_stem + "_noise" + std::to_string(n_idx) , n.size(),
                                    n.dilutionSize(Index::l), n.dilutionSize(Index::s), n.dilutionSize(Index::t), D, traj_);
    STOP_TIMER("IO: read");
    STOP_TIMER("IO: total");
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeRhoComponent(FermionField&        rho_component,
                   DistillationNoise&   n,
                   const uint           n_idx,
                   const uint           D,
                   TimerArray*                                 tarray)   
{
    rho_component = n.makeSource(D, n_idx);
}


template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeDvLapSpinComponent(FermionField&                          component,
                        std::map<Side, uint>                    n_idx,
                        LapPack&                                epack,
                        Side                                    s,
                        uint                                    D,
                        std::map<Side, MDistil::PerambTensor&>  peramb,
                        TimerArray*                                 tarray)
{
    if(isPhi(s))
    {
        if(vectorStem_.at(s).empty())
        {
            makePhiComponent(component , distilNoise_.at(s) , n_idx.at(s) , D , peramb.at(s), epack, tarray);
        }
        else
        {
            loadPhiComponent(component , distilNoise_.at(s) , n_idx.at(s) , D , vectorStem_.at(s), epack, tarray);
        }
    }
    else if(isRho(s))
    {
        makeRhoComponent(component, distilNoise_.at(s) , n_idx.at(s) , D, tarray);
    }
}

// lap-spin blocks have fixed dimensions of (lap-spin dilution size left)x(lap-spin dilution size right)
// and each is identified by the starting posision in time dilution space, (dtL,dtR)
template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeDvLapSpinBlock(std::map<Side, DistilVector&>                      dv,
                               std::map<Side, uint>                     n_idx,
                               LapPack&                                 epack,
                               Side                                     s,
                               uint                                     dt,
                               uint                                     iibatch,
                               TimerArray*                                 tarray,
                               std::map<Side, MDistil::PerambTensor&>   peramb)
{
    uint D_offset = distilNoise_.at(s).dilutionIndex(dt,0,0);    // t is the slowest index
    uint iD_offset = iibatch*dilSizeLS_.at(s);

    for(uint iD=iD_offset ; iD<iD_offset+dilSizeLS_.at(s) ; iD++)
    {
        uint D = (iD - iD_offset) + D_offset ;
        makeDvLapSpinComponent(dv.at(s)[iD], n_idx, epack, s, D, peramb, tarray);
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeDvLapSpinCacheBlock(std::map<Side, DistilVector&>                     dv,
                               uint                                         dv_idx_offset, //offset variable
                               std::map<Side, uint>                         n_idx,
                               LapPack&                                     epack,
                               Side                                         s,
                               std::map<Side, MDistil::PerambTensor&>       peramb,
                               TimerArray*                                 tarray)
{
    for(uint iD=0 ; iD<cacheSize_ ; iD++)
    {
        uint D = iD + dv_idx_offset;
        makeDvLapSpinComponent(dv.at(s)[iD], n_idx, epack, s, D, peramb, tarray);
    }
    // if(isPhi(s) and vectorStem_.at(s).empty())
    // {
    //     LOG(Message) << "*** Breakdown makeDvLapSpinCacheBlock -> makePhiComponent timers:" << std::endl;
    //     LOG(Message) << "* ExtractSliceLocal " << GET_TIMER("makePhiComponent:ExtractSliceLocal")/1e6 << " s"  << std::endl;
    //     LOG(Message) << "* tmp3d_+=evec3d_ " << GET_TIMER("makePhiComponent:ExtractSliceLocal")/1e6 << " s"  << std::endl;
    //     LOG(Message) << "* InsertSliceLocal " << GET_TIMER("makePhiComponent:InsertSliceLocal")/1e6 << " s"  << std::endl;
    // }
    // else if(isRho(s))
    // {
    // }
}

// a batch is composed of several lap-spin blocks
template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeDvLapSpinBatch(std::map<Side, DistilVector&>                      dv,
                               std::map<Side, uint>                     n_idx,
                               LapPack&                                 epack,
                               Side                                     s,
                               std::vector<uint>                        dt_list,
                               std::map<Side, MDistil::PerambTensor&>   peramb,
                               TimerArray*                                 tarray)
{
    for(uint idt=0 ; idt<dt_list.size() ; idt++)
    {
        makeDvLapSpinBlock(dv,n_idx,epack,s,dt_list[idt],idt,tarray,peramb);
    }
    // if(isPhi(s) and vectorStem_.at(s).empty())
    // {
        // LOG(Message) << "*** Breakdown makeDvLapSpinCacheBlock -> makePhiComponent timers:" << std::endl;
        // LOG(Message) << "* ExtractSliceLocal " << GET_TIMER("makePhiComponent:ExtractSliceLocal")/1e6 << " s"  << std::endl;
        // LOG(Message) << "* tmp3d_+=evec3d_ " << GET_TIMER("makePhiComponent:ExtractSliceLocal")/1e6 << " s"  << std::endl;
        // LOG(Message) << "* InsertSliceLocal " << GET_TIMER("makePhiComponent:InsertSliceLocal")/1e6 << " s"  << std::endl;
    // }
    // else if(isRho(s))
    // {
    // }
}


// DV RELATIVE METHODS

// nb: methods making relative components just insert a couple of timeslices (1 for full time dilution) on the field 

// makeRelative methods build distil vectors with reorganised time slices
// (in order to compute multiple time-dilution blocks at different t with a single call of MesonField kernel)
template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeRelativePhiComponent(FermionField&                    phi_component,
                   DistillationNoise&                       n,
                   const uint                               n_idx,
                   const uint                               D,
                   const uint                               delta_t,
                   std::map<Side, MDistil::PerambTensor&>   peramb,
                   Side                                     s,
                   LapPack&                                 epack,
                   TimerArray*                                 tarray)
{
    std::array<uint,3> d_coor = n.dilutionCoordinates(D);
    uint dt = d_coor[Index::t] , dk = d_coor[Index::l] , ds = d_coor[Index::s];

    if(dk==0 and ds==0)
    {
        LOG(Message) << "Phi(T=" << dt << ") : t={" << MDistil::timeslicesDump(n.dilutionPartition(Index::t, dt)) << "}" << std::endl;
    }

    for (uint it : n.dilutionPartition(Index::t, dt))
    {
        //compute at rel_t, insert at t=rel_t (full t-dil only)
        const uint rel_t = (it+delta_t)%nt_; // enabling delta_t!=0 only for non t-dilution at Module
        const uint t     = rel_t;

        if(!vectorStem_.at(s).empty())
        {
            if (delta_t!=0) // would have to change D into something aware of rel_t 
            {
                HADRONS_ERROR(Implementation, "delta_t!=0 not supported when loading vector from disk...");
            }
            loadPhiComponent(tmp4d_, n , n_idx , D , vectorStem_.at(s) , epack, tarray);    // potentially io demanding
        }
        
        const uint nVec = epack.evec.size();
        const uint Nt_first = g_->LocalStarts()[nd_ - 1];
        const uint Nt_local = g_->LocalDimensions()[nd_ - 1];
        if( (rel_t>=Nt_first) and (rel_t<Nt_first+Nt_local) )
        {
            tmp3d_ = Zero();
            for (uint k = 0; k < nVec; k++)
            {
                if(vectorStem_.at(s).empty())
                {
                    std::vector<int> peramb_ts = peramb.at(s).MetaData.timeSources;
                    std::vector<int>::iterator itr_dt = std::find(peramb_ts.begin(), peramb_ts.end(), dt);
                    uint idt_peramb = std::distance(peramb_ts.begin(), itr_dt); //gets correspondent index of dt in the tensor obj 
                    ExtractSliceLocal(evec3d_,epack.evec[k],0,rel_t-Nt_first,nd_ - 1);
                    tmp3d_ += evec3d_ * peramb.at(s).tensor(rel_t, k, dk, n_idx, idt_peramb, ds);
                }
                else
                {
                    ExtractSliceLocal(tmp3d_,tmp4d_,0,rel_t-Nt_first,nd_ - 1); // extracting timeslice from loaded field
                }
            }
            InsertSliceLocal(tmp3d_,phi_component,0,t-Nt_first,nd_ - 1);
        }
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeRelativeRhoComponent(FermionField&    rho_component,
                   DistillationNoise&       n,
                   const uint               n_idx,
                   const uint               D,
                   LapPack&                 epack,
                   TimerArray*                                 tarray)
{
    // LOG(Message) << std::endl << "-- makeRelativeRhoComponent(D= " << D << ", delta_t=" << delta_t << ")" << std::endl;
    std::array<uint,3> d_coor = n.dilutionCoordinates(D);
    uint dt = d_coor[Index::t] , dk = d_coor[Index::l] , ds = d_coor[Index::s];

    typename DistillationNoise::NoiseType   noise(n.getNoise()[n_idx].data(), nt_, epack.eval.size(), Ns);

    if(dk==0 and ds==0)
    {
        LOG(Message) << "Rho(T=" << dt << ") : t={" << MDistil::timeslicesDump(n.dilutionPartition(Index::t, dt)) << "}" << std::endl;
    }
    
    for (uint it : n.dilutionPartition(Index::t, dt))
    {
        //compute at rel_t, insert at t (for relative rho vector, delta_t should have been to 0 already, and so they're the same)
        const uint rel_t = it;
        const uint t     = rel_t;

        // adapt to time dilution!
        const uint Nt_first = g_->LocalStarts()[nd_ - 1];
        const uint Nt_local = g_->LocalDimensions()[nd_ - 1];
        if( (rel_t>=Nt_first) and (rel_t<Nt_first+Nt_local) )
        {
            for (uint ik : n.dilutionPartition(Index::l, dk))
            {
                for (uint is : n.dilutionPartition(Index::s, ds))
                {
                    // LOG(Message) << "block(ik=" << ik << ",is=" << is << ")" << std::endl;
                    
                    START_TIMER("makeRelativeRhoComponent:ExtractSliceLocal");
                    ExtractSliceLocal(evec3d_, epack.evec[ik], 0, rel_t - Nt_first, nd_-1);
                    STOP_TIMER("makeRelativeRhoComponent:ExtractSliceLocal");
                    
                    evec3d_ = evec3d_*noise(it, ik, is);
                    tmp3d_  = Zero();
                    pokeSpin(tmp3d_, evec3d_, is);
                    tmp4d_ = Zero();

                    START_TIMER("makeRelativeRhoComponent:InsertSliceLocal");
                    InsertSliceLocal(tmp3d_, tmp4d_, 0, t - Nt_first, nd_-1);
                    STOP_TIMER("makeRelativeRhoComponent:InsertSliceLocal");

                    START_TIMER("makeRelativeRhoComponent:rho_component+=tmp4d_");
                    rho_component += tmp4d_;
                    STOP_TIMER("makeRelativeRhoComponent:rho_component+=tmp4d_");
                }
            }
        }
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::makeRelativeDvLapSpinBlock(std::map<Side, DistilVector&>              dv,
                               std::vector<uint>                        dt_list,
                               std::map<Side, uint>                     n_idx,
                               LapPack&                                 epack,
                               Side                                     s,
                               const uint                               delta_t,
                               std::map<Side, MDistil::PerambTensor&>   peramb,
                               TimerArray*                                 tarray)
{
    DistilVector& dv_rel = dv.at(s);
    DistillationNoise& noise_rel = distilNoise_.at(s);
    for(uint D_LS=0 ; D_LS<dilSizeLS_.at(s) ; D_LS++)    // reset dv
        dv_rel[D_LS] = Zero();

    for(uint D=0 ; D<noise_rel.dilutionSize() ; D++) //loop over all (dt,dk,ds) compound indices
    {
        std::array<uint,3> d_coor = noise_rel.dilutionCoordinates(D);
        uint dt = d_coor[Index::t] , dk = d_coor[Index::l] , ds = d_coor[Index::s];
        // select the dt's wanted here: they will correspond to different timeslices of a single dt component
        // for a rho field, dt_list should always be the full list, as using subsets will still take the same compute & memory (the whole FermionField needs to be used)
        if( std::count(dt_list.begin(), dt_list.end(), dt)!=0 )
        {
            // will compute Field at D, but insert it at D_LS 
            const uint D_LS = noise_rel.dilutionIndex(0,dk,ds);
            if(isPhi(s))
            {
                makeRelativePhiComponent(dv_rel[D_LS], noise_rel, n_idx.at(s), D, delta_t, peramb, s, epack, tarray);
            }
            else if(isRho(s))
            {
                makeRelativeRhoComponent(dv_rel[D_LS], noise_rel, n_idx.at(s), D, epack, tarray); // delta_t=0 always with rho fields
            }
        }
    }
    // if(isRho(s))
    // {
    //     LOG(Message) << "*** Breakdown makeRelativeDvLapSpinBlock -> makeRelativeRhoComponent timers:" << std::endl;
    //     LOG(Message) << "* ExtractSliceLocal " << GET_TIMER("makeRelativeRhoComponent:ExtractSliceLocal")/1e6 << " s"  << std::endl;
    //     LOG(Message) << "* InsertSliceLocal " << GET_TIMER("makeRelativeRhoComponent:InsertSliceLocal")/1e6 << " s"  << std::endl;
    //     LOG(Message) << "* rho_component+=tmp4d_ " << GET_TIMER("makeRelativeRhoComponent:rho_component+=tmp4d_")/1e6 << " s"  << std::endl;
    // }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::executeRelative(const FilenameFn&                         filenameDmfFn,
                const MetadataFn&                           metadataDmfFn,
                Vector<Tio>&                                bBuf,
                Vector<T>&                                  cBuf,
                std::vector<Gamma::Algebra>                 gamma,
                std::map<Side, DistilVector&>               dv,
                std::map<Side, uint>                        n_idx,
                std::vector<ComplexField>                   ph,
                std::map<Side, std::vector<uint>>           time_dil_source,
                LapPack&                                    epack,
                TimerArray*                                 tarray,
                Side                                        relative_side,
                std::vector<uint>                           delta_t_list,
                std::map<Side, MDistil::PerambTensor&>      peramb)
{
    const uint vol = g_->_gsites;
    //parallel IO info
    const uint i_rank = g_->ThisRank();
    const uint N_ranks = g_->RankCount();
    const uint nExtStr = nExt_*nStr_;
    const uint nExtStrLocal = g_->IsBoss() ? nExtStr/N_ranks 
                                        + nExtStr%N_ranks : nExtStr/N_ranks; // put remainder in boss node

    Side anchored_side = (Side::right==relative_side ? Side::left : Side::right);

    for(uint delta_t : delta_t_list) // this works with delta_t>0 only for full-dilution
    {        
        // ################################################
        // ############ COMPUTE RELATIVE BLOCK ############
        // ################################################
        LOG(Message) << std::endl;
        LOG(Message) << "Computing/loading " << sideToStr.at(relative_side) << " ('relative') distil vector for Tsrc= " << time_dil_source.at(relative_side) << std::endl;
        START_TIMER("distil vectors relative");
        // double start = usecond();
        makeRelativeDvLapSpinBlock(dv, time_dil_source.at(relative_side), n_idx, epack, relative_side, delta_t, peramb, tarray);
        // double stop = usecond();
        // LOG(Message) << "Relative block took " << (double)(stop-start)/1.e6 << " s" << std::endl;
        STOP_TIMER("distil vectors relative");

        // trivial loop unless non-default DISTILVECTOR_TIME_BATCH_SIZE>1
        for (uint ibatchAnchored=0 ; ibatchAnchored<time_dil_source.at(anchored_side).size()/dvBatchSize_ ; ibatchAnchored++)
        { 
            std::vector<uint> batch_dtAnchored;
            batch_dtAnchored = fetchDvBatchIdxs(ibatchAnchored,time_dil_source,anchored_side);
            for (uint idx_dtAnchored=0 ; idx_dtAnchored<batch_dtAnchored.size() ; idx_dtAnchored++)
            {
                uint Tanchored = batch_dtAnchored[idx_dtAnchored];

                LOG(Message) << std::endl;
                if(Side::right==relative_side)
                {
                    LOG(Message) << "------------ " << Tanchored << " X ( " << MDistil::timeslicesDump(time_dil_source.at(relative_side)) << ") ------------" << std::endl; 
                }
                else
                {
                    LOG(Message) << "------------ ( " << MDistil::timeslicesDump(time_dil_source.at(relative_side)) << ") X " << Tanchored << " ------------" << std::endl; 
                }
                LOG(Message) << "Time shift (deltaT) : " << delta_t << std::endl; 
                
                // loop over IO-size blocks within the given time-dilution indices
                LOG(Message) << "- Loop over Lap-spin IO-sized blocks on both sides"  << std::endl; 
                for(uint iRel=0 ; iRel<dilSizeLS_.at(relative_side) ; iRel+=blockSize_)
                {
                    // LOG(Message) << "-- Loop over Lap-spin IO-sized blocks on the " << sideToStr.at(anchored_side) << " ('anchored')"  << std::endl; 
                    for(uint jAnchor=0 ; jAnchor<dilSizeLS_.at(anchored_side) ; jAnchor+=blockSize_)
                    {
                        double flops=.0, bytes=.0, t_kernel=0., nodes=g_->NodeCount(), t_cachecopy=.0;
                        // rel_block_size is the size of the current block (indexed by iRel): is either blockSize_ or the size of the remainder
                        uint rel_block_size = MIN(dilSizeLS_.at(relative_side)-iRel,blockSize_);
                        uint anchor_block_size = MIN(dilSizeLS_.at(anchored_side)-jAnchor,blockSize_);

                        //translate relative/anchored into left/right
                        uint left_block_size = (Side::left==relative_side) ? rel_block_size : anchor_block_size;
                        uint right_block_size = (Side::right==relative_side) ? rel_block_size : anchor_block_size;
                        DistilMatrixSetIo<Tio> block(bBuf.data(), nExtStrLocal , nt_, left_block_size, right_block_size);

                        uint nBlockRel = dilSizeLS_.at(relative_side)/blockSize_ + (((dilSizeLS_.at(relative_side) % blockSize_) != 0) ? 1 : 0);
                        uint nBlockAnchor = dilSizeLS_.at(anchored_side)/blockSize_ + (((dilSizeLS_.at(anchored_side) % blockSize_) != 0) ? 1 : 0);

                        LOG(Message) << std::endl;
                        LOG(Message) << "#### IO-size block " << jAnchor/blockSize_ + nBlockAnchor*iRel/blockSize_ + 1 
                            << "/" << nBlockRel*nBlockAnchor <<  " ####" << std::endl;
                        LOG(Message) << "Lap-spin index ranges [relative, anchor] = ["
                        << iRel << " .. " << iRel+rel_block_size-1 << ", " 
                        << jAnchor << " .. " << jAnchor+anchor_block_size-1 << "]" << std::endl;

                        // loop over cache blocks within the current block
                        LOG(Message) << "-- Loop over Lap-spin cache-sized subblocks on the " << sideToStr.at(anchored_side) << " ('anchored')"  << std::endl; 
                        for(uint jjAnchor=0 ; jjAnchor<anchor_block_size ; jjAnchor+=cacheSize_)
                        {
                            // cached dv computation needs to be done here, so it doesnt repeat unnecessarily on inner loop
                            // ################################################
                            // ############ COMPUTE ANCHORED BLOCK ############
                            // ################################################
                            
                            //offsetting time direction, then the block and the cache coordinates
                            uint dv_idx_anchored_offset = Tanchored*dilSizeLS_.at(anchored_side) + jAnchor+jjAnchor; 
                            
                            START_TIMER("distil vectors cache block (anchored)");
                            // start = usecond();
                            makeDvLapSpinCacheBlock(dv,dv_idx_anchored_offset,n_idx,epack,anchored_side,peramb,tarray);
                            // stop = usecond();
                            STOP_TIMER("distil vectors cache block (anchored)");

                            // LOG(Message) << "-- Breakdown kernel execution for each <relative block X anchored block>:" << std::endl;
                            for(uint iiRel=0 ; iiRel<rel_block_size ; iiRel+=cacheSize_)
                            {                                
                                uint rel_cache_size = MIN(rel_block_size-iiRel,cacheSize_);      
                                uint anchor_cache_size = MIN(anchor_block_size-jjAnchor,cacheSize_);

                                //translate relative/anchored into left/right
                                uint left_cache_size = (Side::left==relative_side) ? rel_cache_size : anchor_cache_size;
                                uint right_cache_size = (Side::right==relative_side) ? rel_cache_size : anchor_cache_size;
                                DistilMatrixSetCache<T> cache(cBuf.data(), nExt_, nStr_, nt_, left_cache_size, right_cache_size);

                                // get offset on relative dv (anchored side is not bc it already contains exactly what will be used)
                                uint left_dv_idx_offset  = (Side::left==relative_side)  ? iRel+iiRel : 0;
                                uint right_dv_idx_offset = (Side::right==relative_side) ? iRel+iiRel : 0;

                                START_TIMER("kernel");
                                double t_kernel_tmp = -usecond();
                                A2Autils<FImpl>::MesonField(cache, &dv.at(Side::left)[left_dv_idx_offset], &dv.at(Side::right)[right_dv_idx_offset], 
                                    gamma, ph, nd_ - 1);
                                t_kernel_tmp += usecond();
                                STOP_TIMER("kernel");
                                
                                uint nCacheRel = rel_block_size/cacheSize_ + (((rel_block_size % cacheSize_) != 0) ? 1 : 0);
                                uint nCacheAnchor = anchor_block_size/cacheSize_ + (((anchor_block_size % cacheSize_) != 0) ? 1 : 0);
                                LOG(Message) << "--- Kernel on subblock "
                                    << iiRel/cacheSize_ + nCacheRel*jjAnchor/cacheSize_ + 1 << "/" << nCacheRel*nCacheAnchor
                                    << " (" << (double)(t_kernel_tmp)/1e6 << "s)" << std::endl;

                                const double flops_tmp = vol*(2*8.0+6.0+8.0*nExt_)*rel_cache_size*anchor_cache_size*nStr_;
                                const double bytes_tmp = vol*(12.0*sizeof(T))*rel_cache_size*anchor_cache_size
                                                        +  vol*(2.0*sizeof(T)*nExt_)*rel_cache_size*anchor_cache_size*nStr_;

                                flops += flops_tmp;
                                bytes += bytes_tmp;

                                //// copy cache to ioblock
                                //translate relative/anchored into left/right
                                uint left_ii   = (Side::left==relative_side)  ? iiRel : jjAnchor;
                                uint right_jj  = (Side::right==relative_side) ? iiRel : jjAnchor;

                                t_cachecopy -= usecond();
                                g_->Barrier();
                                START_TIMER("cache copy");
                                thread_for_collapse(4,iextstr_local,nExtStrLocal,{
                                for(uint t=0;t<nt_;t++)
                                for(uint iii=0;iii<left_cache_size;iii++)
                                for(uint jjj=0;jjj<right_cache_size;jjj++)
                                {
                                    const uint iextstr = iextstr_local + i_rank * nExtStrLocal + (g_->IsBoss() ? 0 : nExtStr%N_ranks );
                                    const uint iext = iextstr/nStr_;
                                    const uint istr = iextstr%nStr_;
                                    block(iextstr_local,t,left_ii+iii,right_jj+jjj) = cache(iext,istr,t,iii,jjj);
                                }
                                });
                                STOP_TIMER("cache copy");
                                g_->Barrier();
                                t_cachecopy += usecond();
                            }
                        }
                        // LOG(Message) << "Kernel perf (flops) " << flops/t_kernel/1.0e3/nodes 
                        //             << " Gflop/s/node " << std::endl;
                        // LOG(Message) << "Kernel perf (read) " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes
                        //             << " GB/s/node "  << std::endl;
                        blockCounter_++;
                        blockFlops_ += flops/t_kernel/1.0e3/nodes ;
                        blockBytes_ += bytes/t_kernel*1.0e6/1024/1024/1024/nodes;

                        LOG(Message) << "IO-size block computation complete" << std::endl;
                        LOG(Message) << "Cache to IO copy complete, waited all tasks (" << (double)t_cachecopy/1e6 << "s)" << std::endl;

                        // io section
                        DistillationNoise& n_rel = distilNoise_.at(relative_side);
                        LOG(Message) << "#### Parallel IO ####" << std::endl;
                        LOG(Message) << "Rank count=" << N_ranks << std::endl;
                        // for(uint it=0 ; it<time_dil_source.at(relative_side).size() ; it++)
                        for (uint dt : time_dil_source.at(relative_side))
                        {
                            for(uint it : n_rel.dilutionPartition(Index::t, dt))
                            {
                                const uint t = (it + delta_t)%nt_;
                                const uint Trelative = dt;
                                
                                const uint blockdimleft = (Side::left==relative_side) ? rel_block_size : anchor_block_size;
                                const uint blockdimright = (Side::right==relative_side) ? rel_block_size : anchor_block_size;
                                DistilMatrixSetTimeSliceIo<Tio> block_relative(bBuf.data(), nExtStrLocal, blockdimleft, blockdimright);
                                std::string dataset_name = std::to_string( (Side::right==relative_side) ? Tanchored : Trelative ) 
                                    + "-" + std::to_string( (Side::right==relative_side) ? Trelative : Tanchored );

                                std::vector<uint> relative_partition = distilNoise_.at(relative_side).dilutionPartition(Index::t,Trelative);
                                std::vector<uint> anchored_partition = distilNoise_.at(anchored_side).dilutionPartition(Index::t,Tanchored);

                                // if either side is rho, check t is in the dil. partition on that side
                                if( !( isRho(relative_side) and std::count(relative_partition.begin(), relative_partition.end(), t)==0  ) 
                                    and !(isRho(anchored_side) and std::count(anchored_partition.begin(), anchored_partition.end(), t)==0) )
                                {
                                    LOG(Message)    << "Saving block " << dataset_name << " , t=" << t << std::endl;
                                    double ioTime = -GET_TIMER("IO: write block");
                                    START_TIMER("IO: total");
    #ifdef HADRONS_DISTIL_PARALLEL_IO
                                    g_->Barrier();
                                    for(uint iextstr_local=0 ; iextstr_local<nExtStrLocal ; iextstr_local++)
                                    {
                                        const uint iextstr = iextstr_local + i_rank * nExtStrLocal + (g_->IsBoss() ? 0 : nExtStr%N_ranks );
                                        const uint iext = iextstr/nStr_;
                                        const uint istr = iextstr%nStr_;

                                        // io object
                                        DistilMatrixIo<HADRONS_DISTIL_IO_TYPE> matrix_io(filenameDmfFn(iext, istr, n_idx.at(Side::left), n_idx.at(Side::right)),
                                                DISTIL_MATRIX_NAME, nt_, dilSizeLS_.at(Side::left), dilSizeLS_.at(Side::right));

                                        //executes once per file
                                        if( ( Tanchored==time_dil_source.at(anchored_side).front() ) and    //first time-dilution idx at one side
                                            ( Trelative==time_dil_source.at(relative_side).front() ) and    // same as above for the other side
                                            ( t==(time_dil_source.at(relative_side).front() + delta_t_list.front())%nt_ ) and    //first time slice
                                            (iRel==0) and (jAnchor==0) )  //first IO block
                                        {
                                            //fetch metadata
                                            DistilMesonFieldMetadata<FImpl> md = metadataDmfFn(iext,istr,n_idx.at(Side::left),n_idx.at(Side::right),
                                                                    fetchDilutionMap(Side::left),fetchDilutionMap(Side::right));
                                            //init file and write metadata
                                            START_TIMER("IO: file creation");
                                            matrix_io.initFile(md, (iextstr==0));
                                            STOP_TIMER("IO: file creation");
                                        }
                                        
                                        //translate relative/anchored into left/right
                                        uint left_i  = (Side::left==relative_side)  ? iRel : jAnchor;
                                        uint right_j  = (Side::right==relative_side)  ? iRel : jAnchor;
                                        START_TIMER("IO: write block");
                                        matrix_io.saveBlock(block_relative, iextstr_local, left_i, right_j, dataset_name, t, blockSize_);
                                        STOP_TIMER("IO: write block");
                                    }
                                    g_->Barrier();
    #else
        HADRONS_ERROR(Implementation, "DistilMesonField serial IO not implemented.");
    #endif              
                                    STOP_TIMER("IO: total");
                                    ioTime    += GET_TIMER("IO: write block");
                                    uint bytesBlockSize  = static_cast<double>(nExt_*nStr_*rel_block_size*anchor_block_size*sizeof(Tio));
                                    double iospeed = bytesBlockSize/ioTime*1.0e6/1024/1024;
                                    LOG(Message)    << "HDF5 IO done " << sizeString(bytesBlockSize) << " in "
                                                    << ioTime  << " us (" << iospeed << " MB/s)" << std::endl;
                                    blockIoSpeed_ += iospeed;
                                }
                            }
                        }
                        LOG(Message) << "IO completed" << std::endl;
                    }
                }
            }
        }
    }
}

template <typename FImpl, typename T, typename Tio>
void DmfComputation<FImpl,T,Tio>
::executeFixed(const FilenameFn&                        filenameDmfFn,
          const MetadataFn&                             metadataDmfFn,
          Vector<Tio>&                                  bBuf,
          Vector<T>&                                    cBuf,
          std::vector<Gamma::Algebra>                   gamma,
          std::map<Side, DistilVector&>                 dv,
          std::map<Side, uint>                          n_idx,
          std::vector<ComplexField>                     ph,
          std::map<Side, std::vector<uint>>             time_dil_source,
          LapPack&                                      epack,
          TimerArray*                                   tarray,
          bool                                          only_diag,
          const uint                                    diag_shift,
          std::map<Side, MDistil::PerambTensor&>        peramb)
{
    const uint vol = g_->_gsites;
    //parallel IO info
    const uint i_rank = g_->ThisRank();
    const uint N_ranks = g_->RankCount();
    const uint nExtStr = nExt_*nStr_;
    const uint nExtStrLocal = g_->IsBoss() ? nExtStr/N_ranks 
                                        + nExtStr%N_ranks : nExtStr/N_ranks; // put remainder in boss node

    // trivial loop unless non-default DISTILVECTOR_TIME_BATCH_SIZE>1
    for (uint ibatchL=0 ; ibatchL<time_dil_source.at(Side::left).size()/dvBatchSize_ ; ibatchL++)
    {
        LOG(Message) << std::endl;
        LOG(Message) << "Computing/loading left distil vector for Tsrc= " << time_dil_source.at(Side::left) << std::endl; 
        std::vector<uint> batch_dtL = fetchDvBatchIdxs(ibatchL,time_dil_source,Side::left);
        START_TIMER("distil vectors left");
        makeDvLapSpinBatch(dv, n_idx, epack, Side::left, batch_dtL, peramb, tarray);
        STOP_TIMER("distil vectors left");

        for (uint idtL=0 ; idtL<batch_dtL.size() ; idtL++)
        {
            uint dtL = batch_dtL[idtL];
            // trivial loop unless non-default DISTILVECTOR_TIME_BATCH_SIZE>1
            for (uint ibatchR=0 ; ibatchR<time_dil_source.at(Side::right).size()/dvBatchSize_ ; ibatchR++)
            {
                std::vector<uint> batch_dtR = fetchDvBatchIdxs(ibatchR,time_dil_source,Side::right,diag_shift);
                for (uint idtR=0 ; idtR<batch_dtR.size() ; idtR++)
                {
                    uint dtR = batch_dtR[idtR];
                    // fetch necessary time slices for this time-dilution block
                    std::map<Side,std::vector<uint>> time_partition = { {Side::left,{}} , {Side::right,{}}};
                    for(auto s : sides)
                    {
                        if(isPhi(s))
                        {
                            time_partition.at(s).resize(nt_);
                            //phi: fill with all time slices so the intersection with a vector X is equal to X
                            std::iota(std::begin(time_partition.at(s)), std::end(time_partition.at(s)), 0); 
                        }
                        else
                        {
                            time_partition.at(s) = distilNoise_.at(s).dilutionPartition(Index::t, s==Side::left ? dtL : dtR);
                        }
                    }
                    std::vector<uint> ts_intersection;
                    std::set_intersection(time_partition.at(Side::left).begin(), time_partition.at(Side::left).end(), 
                                        time_partition.at(Side::right).begin(), time_partition.at(Side::right).end(),
                                        std::back_inserter(ts_intersection));

                    // only execute when partitions have at least one time slice in common; only computes diagonal (up to diag_shift) when onlyDiag is on
                    if( !ts_intersection.empty() and 
                        ( !only_diag or ((dtL+diag_shift)%distilNoise_.at(Side::right).dilutionSize(Index::t)==dtR) ) )
                    {
                        LOG(Message) << std::endl;
                        LOG(Message) << "------------ time-dilution pair " << dtL << "-" << dtR << " ------------------------" << std::endl; 

                        // loop over blocks within the current time-dilution block
                        LOG(Message) << "- Loop over Lap-spin IO-sized blocks on both sides"  << std::endl; 
                        for(uint i=0 ; i<dilSizeLS_.at(Side::left) ; i+=blockSize_) //set according to memory size
                        for(uint j=0 ; j<dilSizeLS_.at(Side::right) ; j+=blockSize_)
                        {
                            double flops=.0, bytes=.0, t_kernel=0., nodes=g_->NodeCount(), t_cachecopy=.0;

                            // iblock_size is the size of the current block (indexed by i); N_i-i is the size of the possible remainder block
                            uint iblock_size = MIN(dilSizeLS_.at(Side::left)-i,blockSize_);
                            uint jblock_size = MIN(dilSizeLS_.at(Side::right)-j,blockSize_);
                            DistilMatrixSetIo<Tio> block(bBuf.data(), nExtStrLocal , nt_, iblock_size, jblock_size);

                            uint nBlocki = dilSizeLS_.at(Side::left)/blockSize_ + (((dilSizeLS_.at(Side::left) % blockSize_) != 0) ? 1 : 0);
                            uint nBlockj = dilSizeLS_.at(Side::right)/blockSize_ + (((dilSizeLS_.at(Side::right) % blockSize_) != 0) ? 1 : 0);

                            LOG(Message) << std::endl;
                            LOG(Message) << "#### IO-size block " 
                            << j/blockSize_ + nBlocki*i/blockSize_ + 1 
                            << "/" << nBlocki*nBlockj << " ####" << std::endl;
                            LOG(Message) << "Lap-spin index ranges [" << i << " .. " 
                            << i+iblock_size-1 << ", " << j << " .. " << j+jblock_size-1 << "]" 
                            << std::endl;

                            LOG(Message) << "-- Loop over Lap-spin cache-sized subblocks on the right"  << std::endl; 
                            // loop over cache blocks within the current block
                            for(uint jj=0 ; jj<jblock_size ; jj+=cacheSize_)
                            {
                                uint dvR_idx_offset = dtR*dilSizeLS_.at(Side::right) +j+jj;
                                START_TIMER("distil vectors right");
                                makeDvLapSpinCacheBlock(dv,dvR_idx_offset,n_idx,epack,Side::right,peramb,tarray);
                                STOP_TIMER("distil vectors right");

                                for(uint ii=0 ; ii<iblock_size ; ii+=cacheSize_)
                                {
                                    uint icache_size = MIN(iblock_size-ii,cacheSize_);      
                                    uint jcache_size = MIN(jblock_size-jj,cacheSize_);

                                    DistilMatrixSetCache<T> cache(cBuf.data(), nExt_, nStr_, nt_, icache_size, jcache_size);
                                    
                                    uint dv_idxL = idtL*dilSizeLS_.at(Side::left) +i+ii;
                                    
                                    START_TIMER("kernel");
                                    double t_kernel_tmp = -usecond();
                                    A2Autils<FImpl>::MesonField(cache, &dv.at(Side::left)[dv_idxL], &dv.at(Side::right)[0], gamma, ph, nd_ - 1); 
                                    t_kernel_tmp += usecond();
                                    STOP_TIMER("kernel");

                                    uint nCacheLeft = iblock_size/cacheSize_ + (((iblock_size % cacheSize_) != 0) ? 1 : 0);
                                    uint nCacheRight = jblock_size/cacheSize_ + (((jblock_size % cacheSize_) != 0) ? 1 : 0);
                                    
                                    LOG(Message) << "--- Kernel on subblock "
                                    << ii/cacheSize_ + nCacheLeft*jj/cacheSize_ + 1 << "/" << nCacheLeft*nCacheRight
                                    << " (" << (double)(t_kernel_tmp)/1e6 << "s)" << std::endl;

                                    const double flops_tmp = vol*(2*8.0+6.0+8.0*nExt_)*icache_size*jcache_size*nStr_;
                                    const double bytes_tmp = vol*(12.0*sizeof(T))*icache_size*jcache_size
                                                        +  vol*(2.0*sizeof(T)*nExt_)*icache_size*jcache_size*nStr_;
                                    
                                    flops += flops_tmp;
                                    bytes += bytes_tmp;
                                    t_kernel += t_kernel_tmp;

                                    // copy cache to ioblock
                                    t_cachecopy -= usecond();
                                    g_->Barrier();
                                    START_TIMER("cache copy");
                                    uint ts_size = ts_intersection.size();
                                    thread_for_collapse(4,iextstr_local,nExtStrLocal,{
                                    for(uint it=0;it<ts_size;it++)
                                    for(uint iii=0;iii<icache_size;iii++)
                                    for(uint jjj=0;jjj<jcache_size;jjj++)
                                    {
                                        const uint iextstr = iextstr_local + i_rank * nExtStrLocal + (g_->IsBoss() ? 0 : nExtStr%N_ranks );
                                        const uint iext = iextstr/nStr_;
                                        const uint istr = iextstr%nStr_;
                                        block(iextstr_local,ts_intersection[it],ii+iii,jj+jjj) = cache(iext,istr,ts_intersection[it],iii,jjj);
                                    }
                                    });
                                    STOP_TIMER("cache copy");
                                    g_->Barrier();
                                    t_cachecopy += usecond();
                                }
                            }
                            // LOG(Message) << "Kernel perf (flops) " << flops/t_kernel/1.0e3/nodes 
                            //             << " Gflop/s/node " << std::endl;
                            // LOG(Message) << "Kernel perf (read) " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes
                            //             << " GB/s/node "  << std::endl;
                            blockCounter_++;
                            blockFlops_ += flops/t_kernel/1.0e3/nodes ;
                            blockBytes_ += bytes/t_kernel*1.0e6/1024/1024/1024/nodes;

                            LOG(Message) << "IO-size block computation complete" << std::endl;
                            LOG(Message) << "Cache to IO copy complete, waited all tasks (" << (double)t_cachecopy/1e6 << "s)" << std::endl;

                            // io section
                            LOG(Message) << "#### Parallel IO ####" << std::endl;
                            LOG(Message) << "Rank count=" << N_ranks << std::endl;
                            for(uint t=0 ; t<nt_ ; t++)
                            {
                                // TODO: generalise to dilution
                                // uint t = ts_intersection[it];
                                std::vector<uint> left_partition = distilNoise_.at(Side::left).dilutionPartition(Index::t,dtL);
                                std::vector<uint> right_partition = distilNoise_.at(Side::right).dilutionPartition(Index::t,dtR);

                                if( !( isRho(Side::left) and std::count(left_partition.begin(), left_partition.end(), t)==0  ) 
                                    and !(isRho(Side::right) and std::count(right_partition.begin(), right_partition.end(), t)==0) )
                                {
                                    // use same buffer but map it differently 
                                    DistilMatrixSetTimeSliceIo<Tio> block_relative(bBuf.data(), nExtStrLocal , iblock_size, jblock_size); 
                                    std::string dataset_name = std::to_string(dtL)+"-"+std::to_string(dtR);
                                    LOG(Message)    << "Saving block block " << dataset_name << " , t=" << t << std::endl;

                                    double ioTime = -GET_TIMER("IO: write block");
                                    START_TIMER("IO: total");
#ifdef HADRONS_DISTIL_PARALLEL_IO
                                    g_->Barrier();
                                    for(uint iextstr_local=0 ; iextstr_local<nExtStrLocal ; iextstr_local++)
                                    {
                                        const uint iextstr = iextstr_local + i_rank * nExtStrLocal + (g_->IsBoss() ? 0 : nExtStr%N_ranks );
                                        const uint iext = iextstr/nStr_;
                                        const uint istr = iextstr%nStr_;

                                        // io object
                                        DistilMatrixIo<HADRONS_DISTIL_IO_TYPE> matrix_io(filenameDmfFn(iext, istr, n_idx.at(Side::left), n_idx.at(Side::right)),
                                                DISTIL_MATRIX_NAME, nt_, dilSizeLS_.at(Side::left), dilSizeLS_.at(Side::right));

                                        //executes once per file
                                        if( ( dtL==time_dil_source.at(Side::left).front() ) and      //first time-dilution idx at one side
                                            ( dtR==((time_dil_source.at(Side::right).front()+diag_shift)%distilNoise_.at(Side::right).dilutionSize(Index::t)) ) and    // same as above for the other side
                                            ( t==ts_intersection.front() ) and    //first time slice
                                            (i==0) and (j==0) )  //first IO block
                                        {
                                            // fetch metadata
                                            DistilMesonFieldMetadata<FImpl> md = metadataDmfFn(iext,istr,n_idx.at(Side::left),n_idx.at(Side::right),
                                                                    fetchDilutionMap(Side::left),fetchDilutionMap(Side::right));
                                            //init file and write metadata
                                            START_TIMER("IO: file creation");
                                            matrix_io.initFile(md, (iextstr==0));
                                            STOP_TIMER("IO: file creation");
                                        }
                                        START_TIMER("IO: write block");
                                        matrix_io.saveBlock(block_relative, iextstr_local, i, j, dataset_name, t, blockSize_);
                                        STOP_TIMER("IO: write block");
                                    }
                                    g_->Barrier();
#else
    HADRONS_ERROR(Implementation, "DistilMesonField serial IO not implemented.");
#endif              
                                    STOP_TIMER("IO: total");
                                    ioTime    += GET_TIMER("IO: write block");
                                    uint bytesBlockSize  = static_cast<double>(nExt_*nStr_*iblock_size*jblock_size*sizeof(Tio));
                                    double iospeed = bytesBlockSize/ioTime*1.0e6/1024/1024;
                                    LOG(Message)    << "HDF5 IO done " << sizeString(bytesBlockSize) << " in "
                                                    << ioTime  << " us (" << iospeed << " MB/s)" << std::endl;
                                    blockIoSpeed_ += iospeed;
                                }
                            }
                            LOG(Message) << "IO completed" << std::endl;
                        }
                    }
                }
            }
        }
    }
}

END_HADRONS_NAMESPACE

#endif // Distil_matrix_hpp_
