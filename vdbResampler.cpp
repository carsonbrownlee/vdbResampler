#include <openvdb/openvdb.h>
#include <openvdb/tools/GridTransformer.h>

#include <stdio.h>
#include <iostream>
#include <sstream>

using GridType = openvdb::FloatGrid;

int main(int argc, char** argv)
{
    if (argc < 5)
    {
        std::cerr << "ERROR.  usage: ./vdb_resampler <filename.vdb> <i axis divisor> <j axis divisor> <k axis divisor>\n";
        return 1;
    }
    std::string filename(argv[1]);  //input vdb file
    float idiv = atof(argv[2]);   //output dimensions = inputdim/(i,j,k)div
    float jdiv = atof(argv[3]);
    float kdiv = atof(argv[4]);

    openvdb::initialize();
    // Create a VDB file object.
    openvdb::io::File file(filename.c_str());
    // Open the file.  This reads the file header, but not any grids.
    file.open();
    // Loop over all grids in the file and retrieve a shared pointer
    // to the one named "LevelSetSphere".  (This can also be done
    // more simply by calling file.readGrid("LevelSetSphere").)
    openvdb::GridBase::Ptr baseGrid;
    for (openvdb::io::File::NameIterator nameIter = file.beginName();
        nameIter != file.endName(); ++nameIter)
    {
        // Read in only the grid we are interested in.
        if (nameIter.gridName() == "density") {
            baseGrid = file.readGrid(nameIter.gridName());
        } else {
            std::cout << "skipping grid " << nameIter.gridName() << std::endl;
        }
    }
    file.close();

    auto bbox = baseGrid->evalActiveVoxelBoundingBox();
    std::cout << "input volume dimensions: " << bbox << std::endl;

    // From the example above, "LevelSetSphere" is known to be a FloatGrid,
    // so cast the generic grid pointer to a FloatGrid pointer.
    openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);

    const int extent[] = {bbox.min().x(), bbox.min().y(), bbox.min().z(),
      bbox.max().x(), bbox.max().y(), bbox.max().z()};
    int iDim = extent[3]-extent[0];
    int jDim = extent[4]-extent[1];
    int kDim = extent[5]-extent[2];
    int maxDim = bbox.maxExtent();
    int maxResample = 32;  //maximum resampled dimension
    int resampled_dimensions[] = {
        int(float(iDim)/idiv),
        int(float(jDim)/jdiv),
        int(float(kDim)/kdiv),
    };
    std::cout << "resampled dimensions: " << resampled_dimensions[0]
      << " " << resampled_dimensions[1]
      << " " << resampled_dimensions[2] << "\n";
    size_t size = resampled_dimensions[0]*resampled_dimensions[1]
      *size_t(resampled_dimensions[2]);
    std::vector<float> resampled_grid(size, 0.f);
    float istart=extent[0];
    float jstart=extent[1];
    float kstart=extent[2];
    float istep=(extent[3]-extent[0])/float(resampled_dimensions[0]-1);
    float jstep=(extent[4]-extent[1])/float(resampled_dimensions[1]-1);
    float kstep=(extent[5]-extent[2])/float(resampled_dimensions[2]-1);
    GridType::ConstAccessor accessor = grid->getConstAccessor();
    float completedCounter = -0.1f;
    for (size_t k=0;k<resampled_dimensions[2]; k++) {
      for (size_t j=0;j<resampled_dimensions[1]; j++) {
        for (size_t i=0;i<resampled_dimensions[0]; i++) {
            const openvdb::Vec3R ijk(istart+float(i)*istep,
               jstart+float(j)*jstep, kstart+float(k)*kstep);
            GridType::ValueType v2 = openvdb::tools::QuadraticSampler::sample(accessor, ijk);
            size_t index = i+j*resampled_dimensions[0]
             + k*size_t(resampled_dimensions[0]*resampled_dimensions[1]);
            resampled_grid[index] = float(v2);
            float compFrac = float(index)/float(size-1);
            if (compFrac >= (completedCounter+0.1f)) {
              std::cout << int(compFrac*100.0f) << "%...\n";
              completedCounter = compFrac;
            }
        }
      }
    }

    std::cout << "writing file...\n";
    std::stringstream fileName;
    fileName << "output_" << resampled_dimensions[0] << "_" << resampled_dimensions[1]
      << "_" << resampled_dimensions[2] << ".raw";
    FILE* oFile;
    oFile = fopen(fileName.str().c_str(), "wb");
    fwrite(resampled_grid.data(), 1, resampled_grid.size()*sizeof(float), oFile);
    fclose(oFile);

    return 0;
 }
