#ifndef DQSDEFORMER_HH
#define DQSDEFORMER_HH

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "dual_quat_cu.hpp"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include <igl/readTGF.h>
#include <igl/directed_edge_parents.h>

using namespace std;
using namespace Tbx;

class DQSdeformer
{
public:
    DQSdeformer();
    void SetBVHname(string fName) {bvhName = fName;}
    void SetTGFname(string fName) {tgfName = fName;}
    void SetWname(string fName) {wName = fName;}
    void Initialize(vector<G4ThreeVector> verts){
        cout<<"========================INITIALIZE DEFORMER======================="<<endl;
        SetVertices(verts);
        cout<<"set "<<U.size()<<" vertices"<<endl;
        ReadTGF(tgfName);
        cout<<"imported "+tgfName+" -> "<<C.rows()<<" centers / "<<BE.rows()<<" bones"<<endl;
        ReadWeight(wName);
        cout<<"imported "+wName+" -> "<<weights.size()<<" weights"<<endl;
        ReadBVH(bvhName);
        cout<<"imported "+bvhName+" -> "<<postureData.size()<<" frames"<<endl;
        cout<<"================================================================="<<endl;
    }
    vector<G4ThreeVector> GetVertices(int FrameNo){
        dual_quat_deformer(postureData[FrameNo]);
        return V;
    }
    G4ThreeVector GetTrans(int FrameNo) {return trans[FrameNo];}
    G4int GetFrameNo() {return postureData.size();}
private:
    void ReadTGF(string fName){
       igl::readTGF(fName,C,BE);
       igl::directed_edge_parents(BE,P);
    }
    void ReadBVH(string fName);
    void ReadWeight(string fName){
        ifstream ifs(fName);
        if(!ifs.is_open()){
            cout<<fName<< " was not open"<<endl;exit(100);
        }
        double tmp, epsilon(1e-10);
        weights.clear();
        for(size_t r=0;r<U.size();r++){
            double wSum(0);
            map<int, double> ww;
            for(int c=0;c<BE.rows();c++){
                ifs>>tmp;
                if(tmp<epsilon) continue;
                ww[c]=tmp;
                wSum += tmp;
            }
            for(auto &w:ww) w.second/=wSum;
            weights.push_back(ww);
        }ifs.close();
    }
    void SetVertices(vector<G4ThreeVector> verts){
        for(size_t i=0;i<verts.size();i++) U.push_back(Point3(verts[i].getX(), verts[i].getY(), verts[i].getZ()));
    }
    void dual_quat_deformer(vector<Dual_quat_cu> dual_quat);
    Transfo transfo_from_eulerYXZ(double ay, double ax, double az )
    {
      Vec3 vx = {1, 0, 0 }, vy = { 0, 1, 0 }, vz = { 0, 0, 1 };
      Mat3 tx = Transfo::rotate(vx,ax*deg).get_mat3();
      Mat3 ty = Transfo::rotate(vy,ay*deg).get_mat3();
      Mat3 tz = Transfo::rotate(vz,az*deg).get_mat3();
      return Transfo(ty*tx*tz);
    }

private:
    vector<Point3> U;
    vector<G4ThreeVector> V;
    Eigen::MatrixXd C;
    Eigen::MatrixXi BE;
    Eigen::VectorXi P;
    vector<vector<Dual_quat_cu>> postureData;
    vector<G4ThreeVector> trans;
    vector<map<int,double>> weights;
    string bvhName, wName, tgfName;
};

#endif // DQSDEFORMER_HH
