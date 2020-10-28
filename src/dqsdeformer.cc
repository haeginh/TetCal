#include "dqsdeformer.hh"

DQSdeformer::DQSdeformer()
{

}

void DQSdeformer::ReadBVH(string fileN){
    vector<vector<Transfo>> data;
    vector<vector<Dual_quat_cu>> conv_data;

    ifstream ifs(fileN);
    if(!ifs.is_open()) {cout<<fileN + " is not open"; exit(200);}

    string dump;
    while(getline(ifs,dump)){
        stringstream ss(dump);
        ss>>dump;
        if(dump=="MOTION") break;
    }
    int frameNo; double frameT;
    ifs>>dump>>frameNo;
    ifs>>dump>>dump>>frameT;
    cout<<"frame #: "<<frameNo<<" / frame time: "<<frameT<<endl;
    data.resize(frameNo);
    double x,y,z;

    trans.clear();
    for(int i=0;i<frameNo;i++){
        data[i].resize(23);
        vector<Transfo> bvhData;
        ifs>>x>>y>>z;
        trans.push_back(G4ThreeVector(x*cm,y*cm,z*cm));

        //collect a row
        for(int j=0;j<19;j++){
            ifs>>y>>x>>z;
            if(j==5) bvhData.push_back(transfo_from_eulerYXZ(y,z,x));
            else if(j==11) bvhData.push_back(transfo_from_eulerYXZ(y,-z,-x));
            else bvhData.push_back(transfo_from_eulerYXZ(y,x,z));
        }
        Transfo tf4 = Transfo::rotate(Vec3(-0.2,-1,0),85*deg);
        Transfo tf10 = Transfo::rotate(Vec3(-0.2,1,0),85*deg);
        Transfo tf6 = Transfo::rotate(Vec3(0,1,0),80*deg);
        Transfo tf12 = Transfo::rotate(Vec3(0,1,0),-80*deg);
        bvhData[4]  = bvhData[4]*tf4;
        bvhData[6]  = bvhData[6]*tf6;
        bvhData[10]  = bvhData[10]*tf10;
        bvhData[12]  = bvhData[12]*tf12;

        Transfo tf0 = bvhData[0];
        for(size_t j=0;j<BE.size();j++){
            int p = BE[j][0];
            if(P[j]<0) {
                data[i][j] = Transfo::identity();
                continue;
            }
            Transfo tf = Transfo::translate(C[p])*bvhData[p]*Transfo::translate(-C[p]);
            data[i][j]=data[i][P[j]]*tf;
        }
        for(Transfo &tf:data[i])
            tf = tf0*tf;
    }

    for(auto d:data){
        vector<Dual_quat_cu> du_quat;
        for(auto tf:d){
            du_quat.push_back(Dual_quat_cu(tf));
        }
        conv_data.push_back(du_quat);
    }
    postureData.clear();
    postureData = conv_data;
}


void DQSdeformer::dual_quat_deformer(vector<Dual_quat_cu> dual_quat)
{
    cout<<"Start deformation.."<<flush;
    V.clear();
    for(int n=0;n<U.size();n++)
    {
        Dual_quat_cu dq_blend;
        bool first(true);
        Quat_cu q0;
        for(auto w:weights[n]){
            if(first){
                dq_blend = dual_quat[w.first] * w.second;
                q0 = dual_quat[w.first].rotation();
                first = false;
                continue;
            }
            if( dual_quat[w.first].rotation().dot( q0 ) < 0.f )
                dq_blend = dq_blend + dual_quat[w.first] * (-w.second);
            else dq_blend = dq_blend + dual_quat[w.first] * w.second;
        }
        // Compute animated position
        Point3 v = dq_blend.transform(U[n]);
        V.push_back(G4ThreeVector(v.x,v.y,v.z));
    }
    cout<<"done"<<endl;
}
