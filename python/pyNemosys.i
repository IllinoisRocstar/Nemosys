%module pyNemosys
%include "std_string.i"
%include "std_vector.i"
%{
#include "meshBase.H"
%}


%template(vectorString) std::vector<std::string>;


class meshBase
{

  // constructors and destructors
  public:

    meshBase():numPoints(0),numCells(0),hasSizeField(0),checkQuality(0);
    ~meshBase();

    static meshBase* Create(std::string fname);
    static std::shared_ptr<meshBase> CreateShared(std::string fname);
    //static std::unique_ptr<meshBase> CreateUnique(std::string fname);
    virtual std::vector<double> getPoint(int id);
    virtual std::map<int, std::vector<double>> getCell(int id);
    virtual std::vector<std::vector<double>> getCellVec(int id);
    vtkSmartPointer<vtkDataSet> getDataSet();
    virtual void setPointDataArray(const char* name,
                                     const std::vector<std::vector<double>>& data);
    virtual void setCellDataArray(const char* name,
                                    const std::vector<std::vector<double>>& data);
    virtual void setCellDataArray(const char* name,
                                    const std::vector<double>& data);
    virtual void unsetPointDataArray(int arrayID);
    virtual void unsetPointDataArray(const char* name);
    virtual void unsetCellDataArray(int arrayID);
    virtual void unsetCellDataArray(const char* name);
    virtual void unsetFieldDataArray(const char* name);
    virtual std::vector<double> getCellLengths();
    virtual std::vector<double> getCellCenter(int cellID);
    vtkSmartPointer<vtkCellLocator> buildLocator();
    virtual void getIntegrationPointsAtCell(int cellID);int transfer(meshBase* target, std::string method,
                 const std::vector<int>& arrayIDs);
    int transfer(meshBase* target, std::string method,
                 const std::vector<std::string>& arrayNames);
    int transfer(meshBase* target, std::string method);
    void generateSizeField(std::string method, int arrayID, double dev_mlt, bool maxIsmin);

    void setSFBool(bool q);
    bool getSFBool();
    int IsArrayName(std::string name);

    void refineMesh(std::string method, int arrayID,
                    double dev_mult, bool maxIsmin, double edge_scale, std::string ofname);
    void refineMesh(std::string method, std::string arrayName,
                    double dev_mult, bool maxIsmin, double edge_scale, std::string ofname);
    void refineMesh(std::string method, double edge_scale, std::string ofname);

    virtual void report();
    int getNumberOfPoints();
    int getNumberOfCells();
    void checkMesh(std::string ofname);

    virtual void write();
    virtual void write(std::string fname);
    void writeMSH(std::ofstream& outputStream);
    void writeMSH(std::string fname);
    void writeMSH(std::ofstream& outputStream, std::string pointOrCell, int arrayID);
    void writeMSH(std::string fname, std::string pointOrCell, int arrayID);
    void writeMSH(std::ofstream& outputStream, std::string pointOrCell, int arrayID,
                    bool onlyVol);
    void writeMSH(std::string fname, std::string pointOrCell, int arrayID,
                    bool onlyVol);

    void setFileName(std::string fname);
    std::string getFileName();
    void setCheckQuality(bool x);
};

int diffMesh(meshBase* mesh1, meshBase* mesh2);
