#include <vector>
#include <map>
#include <cmath>
#include <iostream>
using namespace std;
//valid tag have to be larger than 0

class particle {
private:
    vector<double> crds;
    int tag;

public:
    explicit particle(int t):crds(),tag(t){}
    particle(double x,double y,int t):crds(2),tag(t){crds[0]=x;crds[1]=y;}
    const vector<double>& coords()const{return crds;}
    vector<double>& coords(){return crds;}
    int gettag()const{return tag;}
    bool sametag(const particle& p)const{return tag==p.tag;}
    void set(double x,double y){crds.resize(2);crds[0]=x;crds[1]=y;}
};

class simplex {
private:
    int tag;
    vector<particle*> ptc;
    vector<simplex*> child;
    vector<simplex*> neighbor;
    vector<double> sphere;
    static double tol;
  
public:
    //constructor
    explicit simplex(int t):tag(t),ptc(),child(),neighbor(),sphere(){sphere.reserve(3);}
    simplex(particle&p0,particle&p1,particle&p2,int t):tag(t),ptc(3),child(),neighbor(3,NULL),sphere(3,0.){ptc[0]=&p0;ptc[1]=&p1;ptc[2]=&p2;}

    //set particles
    void set(particle&p0,particle&p1,particle&p2){ptc.resize(3);ptc[0]=&p0;ptc[1]=&p1;ptc[2]=&p2;child.clear();neighbor.resize(3,NULL);sphere.resize(3,0.);}

    //set child
    void addchild(simplex&c){child.push_back(&c);}

    //set neghbor corresponding to edge p0-p1
    void set(const particle&p0,const particle&p1,simplex&s){
        int index[2]={have(p0),have(p1)};
        if(index[0]==-1||index[1]==-1)return;
        neighbor[3-index[0]-index[1]]=&s;
    }
  
    //set the tolerence
    static void set(double t){tol=t;}

    int gettag()const{return tag;}

    const vector<double>& getsphere()const{return sphere;}

    const vector<particle*>& getptc()const{return ptc;}
    vector<particle*>& getptc(){return ptc;}

    //get particle corresponding to edge p0-p1
    particle* getptc(const particle&p0,const particle&p1){
        int index[2]={have(p0),have(p1)};
        if(index[0]==-1||index[1]==-1)return NULL;
        return ptc[3-index[0]-index[1]];
    }

    const vector<simplex*>& getchild()const{return child;}

    //get neighbor corresponding to edge p0-p1
    simplex* getneighbor(const particle&p0,const particle&p1){
        int index[2]={have(p0),have(p1)};
        if(index[0]==-1||index[1]==-1)return NULL;
        return neighbor[3-index[0]-index[1]];
    }

    //utilities
    int have(const particle&p)const{
        for(int i=0;i<ptc.size();i++){
            if(ptc[i]->sametag(p)) return i;
        }
        return -1;
    }

    //1:in sphere  0:on sphere -1:outside sphere
    int insphere(const particle& p);

    //0: in simplex -1:outside simplex 1:on the edge 2:on a point
    int insimplex(const particle& p)const;

    bool sametag(const simplex&tri)const{return tag==tri.tag;}

    void findcenter2D();
};
 

class delaunay {
private:
    vector<simplex> sim;
    vector<particle> ptc;
public:
    explicit delaunay(double large=1e6,int num=100):sim(),ptc(){
        sim.reserve(5*num);
        ptc.reserve(num);
        ptc.push_back(particle(0));
        ptc.push_back(particle(1));
        ptc.push_back(particle(2));
        sim.push_back(simplex(0));
        ptc[0].set(large,0.);
        ptc[1].set(-large,large);
        ptc[2].set(-large,-large);
        sim[0].set(ptc[0],ptc[1],ptc[2]);
    }
    //the pointer to the simplex that contains p; NULL:p on a paritcle
    simplex* locate(const particle& p);
    //tag of new particle
    int add(double x,double y);
    void legalizeedge(simplex& s,particle&p0,particle&p1);
    const vector<particle>& getptc()const{return ptc;}
    const vector<simplex>& getsim()const{return sim;}
    bool validate();
    void merge();
};

