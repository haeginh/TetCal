#ifndef DQSDEFORMER_HH
#define DQSDEFORMER_HH

//#include <Eigen/Geometry>
//#include <Eigen/StdVector>
//#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "dual_quat_cu.hpp"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
//#include <igl/readTGF.h>
//#include <igl/directed_edge_parents.h>

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
        cout<<"imported "+tgfName+" -> "<<C.size()<<" centers / "<<BE.size()<<" bones"<<endl;
        ReadWeight(wName);
        cout<<"imported "+wName+" -> "<<weights.size()<<" weights"<<endl;
        ReadBVH(bvhName);
        cout<<"imported "+bvhName+" -> "<<postureData.size()<<" frames"<<endl;
        cout<<"================================================================="<<endl;
    }
    vector<G4ThreeVector> GetVertices(int FrameNo){
        dual_quat_deformer(postureData[FrameNo]);
//        ofstream ofs("test"+std::to_string(FrameNo)+".pts");
//    ofs<<"#Version 2.0"<<G4endl;
//    ofs<<"Shell type : GENERIC"<<G4endl;
//    ofs<<"Projection Matrix"<<G4endl;
//    ofs<<"1.0000000000 0.0000000000 0.0000000000 0.0000000000"<<G4endl;
//    ofs<<"0.0000000000 1.0000000000 0.0000000000 0.0000000000"<<G4endl;
//    ofs<<"0.0000000000 0.0000000000 1.0000000000 0.0000000000"<<G4endl;
//    ofs<<"0.0000000000 0.0000000000 0.0000000000 1.0000000000"<<G4endl;
//    ofs<<"Number of Vertices : "<<V.size()<<G4endl;
//        for(auto v:V) ofs<<v.getX()<<" "<<v.getY()<<" "<<v.getZ()<<G4endl;
//        ofs.close();
        return V;
    }
    G4ThreeVector GetTrans(int FrameNo) {return trans[FrameNo];}
    G4int GetFrameNo() {return postureData.size();}
private:
    void ReadTGF(string fName){
        ifstream ifs(fName);
        double x, y, z;
        string dump;
        map<int, int> idConv;
        while(ifs>>dump){
            if(dump=="#") break;
            ifs>>x>>y>>z;
            idConv[atoi(dump.c_str())] = C.size();
            C.push_back(Vec3(x*cm,y*cm,z*cm));
//            C.push_back(Vec3(x,y,z));
        }
        int p, d;
        while (ifs>>p>>d){
            BE.push_back({idConv[p],idConv[d]});
        }
        ifs.close();

        for(size_t i=0;i<BE.size();i++){
            bool chk(false);
            for(size_t j=0;j<BE.size();j++){
                if(j==i) continue;
                if(BE[j][1]!=BE[i][0]) continue;
                P.push_back(j);
                chk = true;
                break;
            }
            if(!chk) P.push_back(-1);
        }

//       igl::readTGF(fName,C,BE);
//       igl::directed_edge_parents(BE,P);
//       C = C*cm;
//       cout<<C<<endl<<endl;
//       cout<<BE<<endl<<endl;
//       cout<<P<<endl;
    }
    void ReadBVH(string fName);
    void ReadWeight(string fName){
        ifstream ifs(fName);
        if(!ifs.is_open()){
            cout<<fName<< " was not open"<<endl;exit(100);
        }
        double tmp, epsilon(1e-10);
        weights.clear();
        int count(0);
        for(size_t r=0;r<U.size();r++){
            double wSum(0);
            map<int, double> ww;
            for(size_t c=0;c<BE.size();c++){
                ifs>>tmp;
                if(tmp<epsilon) continue;
                ww[c]=tmp;
                wSum += tmp;
            }
            for(auto &w:ww) w.second/=wSum;
            if(wSum<0.9) count++;
            weights.push_back(ww);
        }ifs.close();
        if(count>0) {cout<<count<<" vertices were not assigned!!"<<endl; getchar();}
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
    vector<Vec3> C;
    vector<vector<int>> BE;
    vector<int> P;
//    Eigen::MatrixXd C;
//    Eigen::MatrixXi BE;
//    Eigen::VectorXi P;
    vector<vector<Dual_quat_cu>> postureData;
    vector<G4ThreeVector> trans;
    vector<map<int,double>> weights;
    string bvhName, wName, tgfName;
};

#endif // DQSDEFORMER_HH
