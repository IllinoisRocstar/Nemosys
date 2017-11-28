#ifndef TRANSFER_H
#define TRANSFER_H

#include <meshBase.H>
#include <ANN.h>
#include <chrono>
/* TODO: This class is responsible for transfering data between meshes
   TODO: Methods will include:
      1) transfer function using several methods for point data and cell data
      2) transfer quality measures

*/
class Timer1 {
private:
  typedef std::chrono::time_point<std::chrono::system_clock> time_t;

public:
  Timer1() : startTime(), stopTime() {}

  time_t start()   { return (startTime = std::chrono::system_clock::now()); }
  time_t stop()    { return (stopTime  = std::chrono::system_clock::now()); }
  double elapsed() { return std::chrono::duration_cast<std::chrono::milliseconds>
                                                      (stopTime-startTime).count(); }

private:
  time_t startTime, stopTime;
};


class Transfer
{
  // constructors and destructors
  public:
    Transfer(){};
    Transfer(meshBase* _source, meshBase* _target):
      source(_source),target(_target){}
    ~Transfer(){}
  
  // transfer methods
  public:
    int run(); 

  private:
    meshBase* source;
    meshBase* target;   

};

int Transfer::run()
{
  if (!(source && target))
  {
    std::cout << "source and target meshes must be initialized" << std::endl;
    exit(1);
  }

  //TODO: cleaning out existing data from target
  //TODO: Implement version for tet meshes using barycentric coords
  //      THIS WILL BE MUCH FASTER

  int nDim, nPnt;
  nDim = 3;
  nPnt = source->getNumberOfPoints();

  // constructing kd-tree from points in source 
  ANNpointArray pntCrd = annAllocPts(nPnt, nDim);
  // filling up vertex coordinate array for the current mesh
  for (int i = 0; i < source->getNumberOfPoints(); ++i)
  {
    std::vector<double> point = source->getPoint(i);
    for (int j = 0; j < nDim; ++j)
      pntCrd[i][j] = point[j];
  } 

  // building kdTree
  ANNkd_tree* kdTree = new ANNkd_tree(pntCrd, nPnt, nDim);

  // transferring point data
  vtkSmartPointer<vtkPointData> pd = vtkSmartPointer<vtkPointData>::New();
  // extracting point data from target
  pd = source->getDataSet()->GetPointData();
  if (pd)
  {
    // finding number of data arrays
    int numArr = pd->GetNumberOfArrays();
    for (int i = 0; i < numArr; ++i)
    { 
      //(pd->GetArrayName(i);
      // extracting data array
      std::vector<std::vector<double>> transferData;
      vtkDataArray* da = pd->GetArray(i);
      // finding dim of data in array 
      int numComponent = da->GetNumberOfComponents();
      transferData.resize(target->getNumberOfPoints());
      

      for (int j = 0; j < target->getNumberOfPoints(); ++j)
      {
        transferData[j].resize(numComponent);
        // getting point from target and setting as query
        std::vector<double> point = target->getPoint(j); 
        ANNpoint queryPt = annAllocPt(nDim);
        for (int k = 0; k < nDim; ++k)
          queryPt[k] = point[k];
        ANNidxArray nnIdx = new ANNidx[1];
        ANNdistArray dists = new ANNdist[1];
        // searching for nearest neighbor of query in source mesh
        kdTree->annkSearch(queryPt,1,nnIdx,dists,0); 
        // if query is a point in source mesh, assign existing data 
        int id = nnIdx[0];
        if (dists[0] == 0)
        {
          double comps[numComponent];
          da->GetTuple(id, comps);
          for (int i = 0; i < numComponent; ++i)
            transferData[j][i] = comps[i];
        }
        else
        {
          double tmp[nDim];
          for (int l = 0; l < nDim; ++l)
            tmp[l] = queryPt[l];
          std::vector<int> commonCells = source->getCellsWithPoint(id);
          for (int k = 0; k < commonCells.size(); ++k)
          {
            int subId;
            double minDist2; // not used
            double pcoords[3];
            int numPointsInCell = source->getDataSet()->GetCell(commonCells[k])
                                  ->GetNumberOfPoints();
            double weights[numPointsInCell];

            int result = source->getDataSet()->GetCell(commonCells[k])->
                         EvaluatePosition(tmp,NULL,subId,pcoords,minDist2,weights);


            if (result == 1)
            {
              for (int h = 0; h < numComponent; ++h)
              {
                for (int m = 0; m < numPointsInCell; ++m)
                {
                  int pntId = source->getDataSet()->GetCell(commonCells[k])->GetPointId(m);
                  double comps[numComponent];
                  da->GetTuple(pntId,comps);
                  transferData[j][h] += comps[h]*weights[m];

                } 
              }
              break;
            }
            else if (result == -1)
            {
              std::cout << "problem encountered evaluating position of point from target"
                        << " mesh with respect to source mesh" << std::endl;
              exit(1);
            }
          }
        }
        delete [] nnIdx;
        delete [] dists;
      }

    //  da->Delete();
      target->setPointDataArray(pd->GetArrayName(i),transferData);
    }  

  }

  delete kdTree;
  annClose();
  
  return 0;

}

#endif
