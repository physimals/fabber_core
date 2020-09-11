/**
 * fabber_rundata_newimage.cc
 *
 * Martin Craig
 *
 * Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "rundata_gifti.h"

#include "easylog.h"
#include "rundata.h"

#include <fslsurface/fslsurfaceio.h>
#include <fslsurface/fslsurface.h>
#include <newmesh/newmesh.h>
#include <newmat.h>

#include <ostream>
#include <string>
#include <vector>

using namespace std;
using NEWMAT::Matrix;
using fslsurface_name::fslSurface;
using fslsurface_name::vertex;
using fslsurface_name::vec3;

FabberRunDataGifti::FabberRunDataGifti(bool compat_options)
    : FabberRunData(compat_options)
{
}

void FabberRunDataGifti::Initialize()
{
    string surface_fname = GetString("surface");
    LOG << "FabberRunDataGifti::Initialize - Reading surface file " << surface_fname << endl;
    fslsurface_name::read_surface(m_surface, surface_fname);
    LOG << "FabberRunDataGifti::Initialize - Vertices=" << m_surface.getNumberOfVertices() << ", faces=" << m_surface.getNumberOfFaces() << endl;
    SetCoordsFromSurface();
}

const Matrix &FabberRunDataGifti::LoadVoxelData(const std::string &filename)
{
    if (m_voxel_data.find(filename) == m_voxel_data.end())
    {
        LOG << "FabberRunDataGifti::LoadVoxelData filename '" + filename << "'" << endl;
        // Load the data file using fslsurface library
        // FIXME should check for presence of file before trying to load.
        //if (!fsl_imageexists(filename))
        //{
        //    throw DataNotFound(filename, "File is invalid or does not exist");
        //}

        fslSurface<float, unsigned int> vertex_data;
        try
        {
            fslsurface_name::read_surface(vertex_data, filename);
        }
        catch (...)
        {
            throw DataNotFound(filename, "Error loading surface data file");
        }

        LOG << "FabberRunDataGifti::LoadVoxelData - Scalars=" << vertex_data.getNumberOfScalarData() 
            << ", vectors=" << vertex_data.getNumberOfVectorData() 
            << endl;

        // We have no way of knowing which array the user wants to use so ensure there is only one
        vector<float> scalars = vertex_data.getScalars(0);
        if (vertex_data.getNumberOfScalarData() != 1) 
        {
            throw DataNotFound(filename, "Expected a single data array from Gifti file - found " + stringify(vertex_data.getNumberOfScalarData()));
        }

        // libfslsurface seems to write timeseries data as scalars with length NV*NT
        // We will use this convention for now but probably should raise this as a enhancement
        // to support timeseries data properly
        unsigned int nv = m_surface.getNumberOfVertices();
        if (scalars.size() % nv != 0) 
        {
            throw DataNotFound(filename, "Data not consistent with surface: " + stringify(scalars.size()) + 
                                         " values for " + stringify(nv) + " vertices");
        }
        unsigned int nt = scalars.size() / nv;
        NEWMAT::Matrix voxel_data(nt, nv);

        // We are making an assumption here about array ordering, but this is a result of 
        // libfslsurface munging a 2D array into 1D
        for (unsigned int v=0; v<nv; v++) 
        {
            for (unsigned int t=0; t<nt; t++) 
            {
                voxel_data(t+1, v+1) = scalars[v*nt+t];
            }
        }

        m_voxel_data[filename] = voxel_data;
    }

    return m_voxel_data[filename];
}

void FabberRunDataGifti::SaveVoxelData(
    const std::string &filename, NEWMAT::Matrix &data, VoxelDataType data_type)
{
    fslSurface<float, unsigned int> output_data;
    unsigned int nv = m_surface.getNumberOfVertices();
    if (data.Ncols() != nv)
    {
        throw FabberInternalError("FabberRunDataGifti::Output data " + filename +
                                  " not consistent with reference surface: " + stringify(data.Ncols()) + 
                                  " values for " + stringify(nv) + " vertices");
    }
    unsigned int nt = data.Nrows();

    vector<float> vertex_data(nt*nv);
    for (unsigned int v=0; v<nv; v++) 
        {
            for (unsigned int t=0; t<nt; t++) 
            {
                vertex_data[v*nt+t] = data(t+1, v+1);
            }
        }

    output_data.addScalars(vertex_data, filename);

    if (filename[0] == '/')
    {
        // Absolute path
        fslsurface_name::writeGIFTI(output_data, filename + ".func.gii");
    }
    else
    {
        // Relative path
        string filepath = GetOutputDir() + "/" + filename;
        fslsurface_name::writeGIFTI(output_data, filepath + ".func.gii");
    }
}

void FabberRunDataGifti::SetCoordsFromSurface()
{
    LOG << "FabberRunDataGifti::Setting coordinates from surface" << endl;

    Matrix coords(3, m_surface.getNumberOfVertices());

    unsigned int idx = 1;
    for (std::vector<vertex<float> >::iterator iter=m_surface.vbegin();
         iter != m_surface.vend(); ++iter) 
    {
        coords(1, idx) = iter->x;
        coords(2, idx) = iter->y;
        coords(3, idx) = iter->z;
        ++idx;
    }

    SetVoxelCoords(coords);
}

static float dot(const vec3<float> &v1, const vec3<float> &v2)
{
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

static vec3<float> cross(const vec3<float> &v1, const vec3<float> &v2)
{
    return vec3<float>(
        v1.y * v2.z - v1.z * v2.y, 
        v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x
    );
}

void FabberRunDataGifti::GetNeighbours(std::vector<std::vector<int> > &neighbours, 
                                       std::vector<std::vector<int> > &neighbours2,
                                       std::vector<std::vector<double> > &weightings)
{
    LOG << "FabberRunDataGifti::Getting nearest neigbours and second neighbours" << endl;

    unsigned int nv = m_surface.getNumberOfVertices();
    unsigned int nt = m_surface.getNumberOfFaces();
    neighbours.clear();
    neighbours2.clear();
    neighbours.resize(nv);
    neighbours2.resize(nv);

    // Iterate over faces
    for (unsigned int t=0; t<nt; t++)
    {
        vec3<unsigned int> trig = m_surface.getFace(t, 3);

        // Iterate over vertices in a triangle - each vertex is connected to each other vertex
        // so add neighbours for all combinations unless they are already in the neighbours list
        // (no duplicates in the nearest neighbour lists)
        for (int v1=0; v1<3; v1++) 
        {
            if (std::find(neighbours[trig.x].begin(), neighbours[trig.x].end(), trig.y + 1) == neighbours[trig.x].end())
                neighbours[trig.x].push_back(trig.y + 1);
            if (std::find(neighbours[trig.x].begin(), neighbours[trig.x].end(), trig.z + 1) == neighbours[trig.x].end())
                neighbours[trig.x].push_back(trig.z + 1);
            if (std::find(neighbours[trig.y].begin(), neighbours[trig.y].end(), trig.x + 1) == neighbours[trig.y].end())
                neighbours[trig.y].push_back(trig.x + 1);
            if (std::find(neighbours[trig.y].begin(), neighbours[trig.y].end(), trig.z + 1) == neighbours[trig.y].end())
                neighbours[trig.y].push_back(trig.z + 1);
            if (std::find(neighbours[trig.z].begin(), neighbours[trig.z].end(), trig.x + 1) == neighbours[trig.z].end())
                neighbours[trig.z].push_back(trig.x + 1);
            if (std::find(neighbours[trig.z].begin(), neighbours[trig.z].end(), trig.y + 1) == neighbours[trig.z].end())
                neighbours[trig.z].push_back(trig.y + 1);
        }
    }

    // Search for Neighbours-of-neighbours, excluding self,
    // but including duplicates if there are two routes to get there (diagonally connected)
    // Note that voxel indices start at 1 because of NEWMAT
    for (unsigned int vid=1; vid <= nv; vid++)
    {   
        //cout << "Voxel " << (vid-1);
        // Go through the list of neighbours for each voxel.
        for (unsigned int n1 = 0; n1 < neighbours[vid-1].size(); n1++)
        {
            //cout << " " << (neighbours[vid-1][n1] - 1);
            // n1id is the voxel index (starting at 1) of the neighbour
            unsigned int n1id = neighbours[vid-1][n1];
            int found_self_as_second_neighbour = 0;
            // Go through each of it's neighbours. Add each, apart from original voxel
            for (unsigned int n2=0; n2<neighbours[n1id-1].size(); n2++)
            {
                unsigned int n2id = neighbours[n1id-1][n2];
                if (n2id != vid)
                {
                    neighbours2[vid-1].push_back(n2id);
                }
                else
                {
                    found_self_as_second_neighbour++;
                }
            }

            if (found_self_as_second_neighbour != 1)
            {
                throw FabberInternalError("Each voxel's neighbour must have "
                                          "the original voxel as a neighbour exactly once");
            }
        }
        //cout << endl;
    }

    weightings.clear();

    // Hold the Voroni masses for each vertex
    std::vector<float> vertex_masses(nv, 0);

    // Initialize the weightings list for each vertex with a vector of zeros for each neighbour
    // plus an extra entry for the vertex self-weight
    for (unsigned int vidx=0; vidx<nv; vidx++)
    {
        int nn = neighbours[vidx].size();
        weightings.push_back(vector<double>(nn+1, 0));
    }

    // This code calculates the weighting for each neighbour of a vertex
    // using the mesh Laplacian
    // Python code to duplicate:
    // vals = []
    // rows = []
    // cols = []
    // avals = []
    // arows = []
    // acols = []
    // for idx, t in enumerate(f):
    for (unsigned int tidx=0; tidx<nt; tidx++)
    {
        vec3<unsigned int> trig = m_surface.getFace(tidx, 3);
        vector<unsigned int> t;
        t.push_back(trig.x);
        t.push_back(trig.y);
        t.push_back(trig.z);

        //     # Triangle vertices (real co-ordinates)
        //     tv = np.array([
        //         v[t[0]],
        //         v[t[1]],
        //         v[t[2]],
        //     ])
        vector<vec3<float> > tv;
        tv.push_back(m_surface.getVertexCoord(t[0]));
        tv.push_back(m_surface.getVertexCoord(t[1]));
        tv.push_back(m_surface.getVertexCoord(t[2]));

        //     # Triangle edges - edge 0 is opposite vertex 0
        //     # with consistent orientation
        //     te = np.array([
        //         tv[2] - tv[1],
        //         tv[0] - tv[2],
        //         tv[1] - tv[0],
        //     ])
        vector<vec3<float> > te;
        te.push_back(tv[2] - tv[1]);
        te.push_back(tv[0] - tv[2]);
        te.push_back(tv[1] - tv[0]);
            
        //     # Dot products of edges at each vertex (negative
        //     # to capture interior angle of triangle)
        //     td = np.array([
        //         -np.dot(te[2], te[1]),
        //         -np.dot(te[0], te[2]),
        //         -np.dot(te[1], te[0]),
        //     ])
        vector<float> td;
        td.push_back(-dot(te[2], te[1]));
        td.push_back(-dot(te[0], te[2]));
        td.push_back(-dot(te[1], te[0]));

        //     cross = np.linalg.norm(np.cross(te[0], te[1]))
        //     ta = cross / 2
        //     obtuse = td[0]*td[1]*td[2] < 0
        float cross_norm = cross(te[0], te[1]).norm();
        float ta = cross_norm/2;
        bool obtuse = td[0]*td[1]*td[2] < 0;

        //     for e in range(3):
        for (unsigned int e=0; e<3; e++) 
        {
            // Consider edge of triangle opposite vertex 'e'
            // cot of angle at vertex 'e' is dot product of other
            // two edges divided by cross product
        //         cot = td[e] / cross
        // vi1, vi2, vi3 = t[e], t[(e+1) % 3], t[(e+2) % 3]
            int vi2 = t[(e+1) % 3];
            int vi3 = t[(e+2) % 3];
            float cot = td[e] / cross_norm;
            
        //         if obtuse:
        //             # Obtuse triangles use barycentric area
        //             # (as in http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.37.8894&rep=rep1&type=pdf)
        //             #avals.append(ta / 3)
        //             if td[e] < 0:
        //                 # Vertex with the obtuse angle gets half the area
        //                 # This is as done in IGL (see http://www.alecjacobson.com/weblog/?p=874)
        //                 avals.append(ta / 2)
        //             else:
        //                 avals.append(ta / 4)
        //             arows.append(t[e])
        //             acols.append(t[e])
        //         else:
        //             # For non obtuse triangles also add the Voronoi area
        //             esqlen = np.sum(np.square(te[e]))
        //             aval = cot * esqlen / 8
        //             avals.append(aval)
        //             avals.append(aval)
        //             arows.append(vi2)
        //             acols.append(vi2)
        //             arows.append(vi3)
        //             acols.append(vi3)
            // Calculate the mass contribution for this triangle. If not obtuse we use the
            // Voroni area, otherwise we give half the area to the vertex with the obtuse
            // angle and quarter to the other two vertices
            // (this is as done in IGL - see http://www.alecjacobson.com/weblog/?p=874)
            if (obtuse)
            {
                // For obtuse triangle contribution is to the vertex
                // opposite the edge depending on whether this vertex
                // has the obtuse angle or not
                if (cot < 0)
                {
                    vertex_masses[t[e]] += ta / 2;
                }
                else
                {
                    vertex_masses[t[e]] += ta / 4;
                }
            }
            else
            {
                // Voronoi area is associated with this edge
                // is applied to the two end vertices of the edge
                float sqnorm = te[e].norm();
                float voronoi_area = cot * sqnorm * sqnorm / 8;
                vertex_masses[vi2] += voronoi_area;
                vertex_masses[vi3] += voronoi_area;
            }

            // Calculate the contribution to the weightings (equivalent to the cotangent
            // matrix). Each edge contributes in 4 places:
            //   (v2, v2), (v3, v3), (v2, v3), (v3, v2)
            // where v2, v3 are the vertices forming the edge. We store the off-diagonal
            // weights in a list for each voxel corresponding to the neighbour voxels
            // in the neighbours list. We add on to the end the diagonal entry for
            // the voxel

        //         vals.append(0.5 * cot)
        //         vals.append(0.5 * cot)
        //         vals.append(-0.5 * cot)
        //         vals.append(-0.5 * cot)
        //         rows.append(vi2)
        //         rows.append(vi3)
        //         rows.append(vi2)
        //         rows.append(vi3)
        //         cols.append(vi2)
        //         cols.append(vi3)
        //         cols.append(vi3)
        //         cols.append(vi2)

            vector<int> v_neighbours = neighbours[vi2];
            for (unsigned int i=0; i<v_neighbours.size(); i++)
            {
                if (v_neighbours[i] == vi3+1) weightings[vi2][i] -= 0.5*cot;
            }
            weightings[vi2].back() += 0.5 * cot;
            v_neighbours = neighbours[vi3];
            for (unsigned int i=0; i<v_neighbours.size(); i++) 
            {
                if (v_neighbours[i] == vi2+1) weightings[vi3][i] -= 0.5*cot;
            }
            weightings[vi3].back() += 0.5 * cot;
        }
    }

    // Finally need to scale weightings by corresponding vertex masses    
    for (unsigned int vidx=0; vidx<nv; vidx++)
    {
        for (unsigned int nidx=0; nidx<weightings[vidx].size(); nidx++)
        {
            weightings[vidx][nidx] /= vertex_masses[vidx];
        }
    }
}

/*
void FabberRunDataGifti::GetWeightings2(std::vector<std::vector<double> > &weightings,
                                       std::vector<std::vector<int> > &neighbours,
                                       const NEWMAT::Matrix &m_laplacian)
{
    LOG << "FabberRunDataGifti::Getting Laplacian weightings for neighbours" << endl;

    int nv = m_surface.getNumberOfVertices();
    weightings.clear();
    weightings.resize(nv);

    // iterate over each vertex
    for (unsigned int vid=1; vid <= nv; vid++)
    {   
        ColumnVector v_weights = m_laplacian.Column(vid);
        // iterate over neighbours of this vertex
        // this only assigns weightings for neighbouring values, not own weighting
        for (unsigned int n2=0; n2<neighbours[vid-1].size(); n2++)
        {
            // get weightings for vertex/neighbour combination from Laplacian matrix
            double weight = v_weights.element(n2);
            weightings[vid-1].push_back(weight);
        }
        // also add the weighting for this vertex
        double weight = v_weights.element(vid-1);
        weightings[vid-1].push_back(weight);
    }
}*/