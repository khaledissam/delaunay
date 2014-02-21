#include "delaunay.hpp"
double simplex::tol=0.0;


void
simplex::findcenter2D()
{
    double x[3],y[3],r[3];
    for(int i=0;i<3;i++)
    {
        x[i]=ptc[i]->coords()[0];
        y[i]=ptc[i]->coords()[1];
        r[i]=x[i]*x[i]+y[i]*y[i];
    }
    double det = 2.*((x[0]-x[1])*(y[0]-y[2])-(x[0]-x[2])*(y[0]-y[1]));
    if(abs(det)<=tol) {
        sphere[2]=-1.;
        return;
    }
    sphere[0] = ((y[0]-y[2])*(r[0]-r[1])-(y[0]-y[1])*(r[0]-r[2]))/det;
    sphere[1] = ((x[0]-x[1])*(r[0]-r[2])-(x[0]-x[2])*(r[0]-r[1]))/det;
    sphere[2] = (sphere[0]-x[0])*(sphere[0]-x[0])+(sphere[1]-y[0])*(sphere[1]-y[0]);
}

int 
simplex::insphere(const particle& p)
{
    //1:in sphere  0:on sphere -1:outside sphere
    const double& xp=p.coords()[0]; const double& yp=p.coords()[1];
    findcenter2D();
    if(sphere[2]==-1) return 1;
  
    double in = sphere[2]-(sphere[0]-xp)*(sphere[0]-xp)-(sphere[1]-yp)*(sphere[1]-yp);
    in/=sphere[2];
    if (in>tol) 
        return 1;
    else if (in<-tol) 
        return -1;
    else 
        return 0;
}

int
simplex::insimplex(const particle& p) const
{
    //0: in simplex -1:outside simplex 1:on an edge 2:on a particle
    if(sphere[2]==-1) return -1;
    const double& x=p.coords()[0]; const double& y=p.coords()[1];
    const double& x1=ptc[0]->coords()[0]; const double& y1=ptc[0]->coords()[1];
    const double& x2=ptc[1]->coords()[0]; const double& y2=ptc[1]->coords()[1];
    const double& x3=ptc[2]->coords()[0]; const double& y3=ptc[2]->coords()[1];
  
    double det = (x1-x3)*(y2-y3)-(x2-x3)*(y1-y3);
    double lambda[3];
    lambda[0] = ((x-x3)*(y2-y3)-(x2-x3)*(y-y3))/det;
    lambda[1] = ((x1-x3)*(y-y3)-(x-x3)*(y1-y3))/det;
    lambda[2] = 1.-lambda[0]-lambda[1];
  
    int res = 0;
    for(int i=0;i<3;i++)
    {
        if(lambda[i]<-tol) return -1;
        if(abs(lambda[i])<=tol)res++;
    }
    return res;
}

simplex*
delaunay::locate(const particle& p)
{
    //the pointer to the simplex that contains p; NULL:p on a paritcle
    simplex* sptr=&sim[0];
    while(true){
        const vector<simplex*>& child=sptr->getchild();
        if(child.empty()) {return sptr;}
        for(int i=0;i<child.size();i++){
            int res = child[i]->insimplex(p);
            if(res==2) return NULL;
            if(res!=-1){
                sptr=child[i];
                break;
            }
        }
    }
}


int
delaunay::add(double x,double y)
{
    //return the tag of the new particle
    int newptc=ptc.size();
    ptc.push_back(particle(newptc));
    ptc.back().set(x,y);
    cout<<"Add new particle "<<ptc.back().gettag()<<endl;

    simplex*sptr=locate(ptc.back());
    if(sptr==NULL) {
        //Coincide with an existing particle!Discard it!
        ptc.pop_back();
        return -1;
    }
    //split the sptr into three parts
    int newsim[3];
    newsim[0]=sim.size();
    newsim[1]=newsim[0]+1;
    newsim[2]=newsim[1]+1;
    sim.push_back(simplex(newsim[0]));
    sim.push_back(simplex(newsim[1]));
    sim.push_back(simplex(newsim[2]));
    vector<particle*> ps=sptr->getptc();

    for(int i=0;i<3;i++){
        int j=i+1,k=i+2;
        if(j>2)j-=3;if(k>2)k-=3;
        sptr->addchild(sim[newsim[i]]);
        sim[newsim[i]].set(*ps[j],*ps[k],ptc.back());
        sim[newsim[i]].set(ptc.back(),*ps[k],sim[newsim[j]]);
        sim[newsim[i]].set(ptc.back(),*ps[j],sim[newsim[k]]);
        simplex* neig=sptr->getneighbor(*ps[j],*ps[k]);
        if(neig!=NULL){
            sim[newsim[i]].set(*ps[j],*ps[k],*neig);
            neig->set(*ps[j],*ps[k],sim[newsim[i]]);
        }
    }
    cout<<"Dividing the triangle "<<sptr->gettag()<<" into three triangles: \n";
    for(int i=0;i<3;i++){
        cout<<sim[newsim[i]].gettag()<<" with";
        for(int j=0;j<3;j++){
            cout<<" "<<sim[newsim[i]].getptc()[j]->gettag();
        }
        cout<<endl;
    }

    //legaize three edges
    for(int i=0;i<3;i++){
        int j=i+1,k=i+2;
        if(j>2)j-=3;if(k>2)k-=3;
        legalizeedge(sim[newsim[i]],*ps[j],*ps[k]);
    }
}

void
delaunay::legalizeedge(simplex& s,particle&p0,particle&p1)
{
    cout<<"Legalizing "<<s.gettag()<<" with edge ";
    cout<<p0.gettag()<<"-"<<p1.gettag()<<endl;

    simplex* neig=s.getneighbor(p0,p1);
    if(neig==NULL) {
        cout<<"triangle "<<s.gettag()<<" doesn't have neighbor "<<endl;
        return;
    }

    particle* neigp=neig->getptc(p0,p1);
    if(s.insphere(*neigp)<=0) {
        cout<<"triangle "<<s.gettag()<<" doesn't contain particle "<<neigp->gettag()<<endl;
        return;
    }

    int newsim[2]={sim.size(),sim.size()+1};
    particle*p[2]={&p0,&p1};
    particle* ownp=s.getptc(p0,p1);
    sim.push_back(simplex(newsim[0]));
    sim.push_back(simplex(newsim[1]));
    cout<<"triangle "<<s.gettag()<<" and "<<neig->gettag();
    cout<<" are divided into triangle "<<sim[newsim[0]].gettag()<<" and "<<sim[newsim[1]].gettag()<<endl;
    for(int i=0;i<2;i++){
        int j=i+1;if(j>1)j-=2;
        neig->addchild(sim[newsim[i]]);
        s.addchild(sim[newsim[i]]);
        sim[newsim[i]].set(*p[i],*ownp,*neigp);
        sim[newsim[i]].set(*ownp,*neigp,sim[newsim[j]]);
    
        simplex* neig2=s.getneighbor(*p[i],*ownp);
    
        if(neig2!=NULL){
            sim[newsim[i]].set(*ownp,*p[i],*neig2);
            neig2->set(*ownp,*p[i],sim[newsim[i]]);
        }
        neig2=neig->getneighbor(*p[i],*neigp);
        if(neig2!=NULL){
            sim[newsim[i]].set(*neigp,*p[i],*neig2);
            neig2->set(*neigp,*p[i],sim[newsim[i]]);
        }
    }
    for(int i=0;i<2;i++){
        legalizeedge(sim[newsim[i]],*p[i],*neigp);
    }
}

bool
delaunay::validate()
{
    for(int i=0;i<sim.size();i++){
        if(!sim[i].getchild().empty()) continue;
        for(int j=0;j<ptc.size();j++){
            if(sim[i].have(ptc[j])!=-1)continue;
            if(sim[i].insphere(ptc[j])>0){
                cout<<"triangle "<<sim[i].gettag()<<" contains particle "<<ptc[j].gettag()<<endl;
                return false;
            }
        }
    }
    return true;
}

