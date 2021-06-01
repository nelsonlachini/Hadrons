/*
 * NamedTensorNoMetadata.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 *  Author: Felix Erben <ferben@ed.ac.uk>
 *  Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: ferben <ferben@debian.felix.com>
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

#ifndef Hadrons_NamedTensorNoMetadata_hpp_
#define Hadrons_NamedTensorNoMetadata_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 NamedTensorNoMetadata contains:
 1) Name of the tensor. By default, this is the tag name used for save / load
 2) Eigen::Tensor of type Scalar_ and rank NumIndices_ (row-major order)
 3) Name for each index
 They can be persisted to / restored from disk. During restore, these validations are performed:
   1) Tensor dimensionality must match
   2) IndexNames are validated against current values
   3) If the tensor has non-zero size, the tensor being loaded must have same extent in each dimension
 ******************************************************************************/

extern const std::string NamedTensorNoMetadataFileExtension;

template<typename Scalar_, int NumIndices_>
class NamedTensorNoMetadata : Serializable
{
public:
    using Scalar = Scalar_;
    static constexpr int NumIndices = NumIndices_;
    using ET = Eigen::Tensor<Scalar_, NumIndices_, Eigen::RowMajor>;
    using Index = typename ET::Index;
    using Traits = Grid::EigenIO::Traits<ET>;
    protected:
      Grid::Vector<typename Traits::scalar_type> deviceBuf;
    public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(NamedTensorNoMetadata,
                                    Eigen::TensorMap<ET>,     tensor,
                                    std::vector<std::string>, IndexNames );

    // Name of the object and Index names as set in the constructor
    const std::string                          &Name_;
    const std::array<std::string, NumIndices_> &DefaultIndexNames_;

    // Return the product of a variable number of dimensions
    Eigen::Index NamedTensorNoMetadataSize() { return 1; }
    template<typename... IndexTypes>
    inline Eigen::Index NamedTensorNoMetadataSize(Eigen::Index firstDimension, IndexTypes... otherDimensions)
    {
        if( sizeof...(otherDimensions) )
            firstDimension *= NamedTensorNoMetadataSize( otherDimensions... );
        return firstDimension;
    }

    // Construct a named tensor explicitly specifying size of each dimension
    template<typename... IndexTypes>
    NamedTensorNoMetadata(const std::string &Name,
                                                      const std::array<std::string, NumIndices_> &indexNames,
                                                      Eigen::Index firstDimension, IndexTypes... otherDimensions)
    : deviceBuf( Traits::count * NamedTensorNoMetadataSize(firstDimension, otherDimensions...) ),
      tensor(reinterpret_cast<Scalar *>(&deviceBuf[0]), firstDimension, otherDimensions...),
      IndexNames{indexNames.begin(), indexNames.end()}, Name_{Name}, DefaultIndexNames_{indexNames}
    {
        if(sizeof...(otherDimensions) + 1 != NumIndices_)
	{
            HADRONS_ERROR(Argument, "NamedTensorNoMetadata: dimensions != tensor rank");
	}
    }

    // Do my index names match the default for my type?
    template<typename array_or_vector_of_string>
    bool ValidateIndexNames( const array_or_vector_of_string &CheckNames ) const
    {
        return IndexNames.size() == CheckNames.size() && std::equal( IndexNames.begin(), IndexNames.end(), CheckNames.begin(),
            [](const std::string &s1, const std::string &s2)
            {
                 return s1.size() == s2.size() && std::equal( s1.begin(), s1.end(), s2.begin(),
                     [](const char & c1, const char & c2)
                     { return c1 == c2 || std::toupper(c1) == std::toupper(c2); }); // case insensitive
            });
    }
    bool ValidateIndexNames() const { return ValidateIndexNames(DefaultIndexNames_); }

#ifdef HAVE_HDF5
    using Default_Reader = Grid::Hdf5Reader;
    using Default_Writer = Grid::Hdf5Writer;
#else
    using Default_Reader = Grid::BinaryReader;
    using Default_Writer = Grid::BinaryWriter;
#endif
    
    void write(const std::string &FileName, const std::string &Tag) const
    {
        std::string FileName_{FileName};
        FileName_.append( NamedTensorNoMetadataFileExtension );
        LOG(Message) << "Writing " << Name_ << " to file " << FileName_ << " tag " << Tag << std::endl;
        Default_Writer w( FileName_ );
        write( w, Tag, *this );
    }
    void write(const std::string &FileName) const { return write(FileName, Name_); }

    // Read tensor.
    // Validate:
    //  1) index names (if requested)
    //  2) index dimensions (if they are non-zero when called)
    template<typename Reader> void read(Reader &r, bool bValidate, const std::string &Tag)
    {
        // Grab index names and dimensions
        std::vector<std::string> OldIndexNames{std::move(IndexNames)};
        const typename ET::Dimensions OldDimensions{tensor.dimensions()};
        read(r, Tag, *this);
        const typename ET::Dimensions & NewDimensions{tensor.dimensions()};
        for (int i = 0; i < NumIndices_; i++)
            if(OldDimensions[i] && OldDimensions[i] != NewDimensions[i])
            {
                HADRONS_ERROR(Size,"NamedTensorNoMetadata::read dimension size");
            }
        if (bValidate && !ValidateIndexNames(OldIndexNames))
        {
            HADRONS_ERROR(Definition,"NamedTensorNoMetadata::read dimension name");
        }
    }
    template<typename Reader> void read(Reader &r, bool bValidate = true) { read(r, bValidate, Name_); }

    inline void read (const std::string &FileName, bool bValidate, const std::string &Tag)
    {
        Default_Reader r(FileName + NamedTensorNoMetadataFileExtension);
        read(r, bValidate, Tag);
    }
    inline void read (const std::string &FileName, bool bValidate= true) { return read(FileName, bValidate, Name_); }
};

/******************************************************************************
 Common elements for distillation
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MDistil)

//Eigenvectors of the Laplacian
using LapEvecs = Grid::Hadrons::EigenPack<LatticeColourVector>;

class PerambTensorNoMetadata : public NamedTensorNoMetadata<SpinVector, 6>
{
    public:
    static const std::string                Name__;
    static const std::array<std::string, 6> DefaultIndexNames__;
    // Construct a named tensor explicitly specifying size of each dimension
    template<typename... IndexTypes>
    PerambTensorNoMetadata(Eigen::Index nT, Eigen::Index nVec, Eigen::Index LI, Eigen::Index nNoise, Eigen::Index nT_inv, Eigen::Index SI)
    : NamedTensorNoMetadata{Name__, DefaultIndexNames__, nT, nVec, LI, nNoise, nT_inv, SI} {}
};

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_NamedTensorNoMetadata_hpp_
