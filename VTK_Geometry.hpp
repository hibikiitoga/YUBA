#include <iostream>
#include <vector>
#include "structure.hpp"
#include "Vector3Dedit.hpp"

//CELL info.
std::array<double, 3> Cell_length;
Vector3D center_of_cell;
bool cell_f=false;
bool scalar_f=false;
int Devide_N=2;
int Devide_tube=50;

constexpr double D_ERROR=0.000001;

//make Regular_Icosahedron
void RI(std::vector<Coordinate>& vertex_list,std::vector<Triangle>& triangle_list)
{
   #define Fg 1
   vertex_list.clear();triangle_list.clear();
   constexpr double G=(1+2.2360679774997896)/2.0;
   Coordinate tmp;
   n=0;
   tmp.v.x= Fg;   tmp.v.y=  G;   tmp.v.z=  0;tmp.n=n++;     vertex_list.push_back(tmp);//0
   tmp.v.x=-Fg;   tmp.v.y=  G;   tmp.v.z=  0;tmp.n=n++;     vertex_list.push_back(tmp);//1
   tmp.v.x=-Fg;   tmp.v.y= -G;   tmp.v.z=  0;tmp.n=n++;     vertex_list.push_back(tmp);//2
   tmp.v.x= Fg;   tmp.v.y= -G;   tmp.v.z=  0;tmp.n=n++;     vertex_list.push_back(tmp);//3
   tmp.v.x=  0;   tmp.v.y= Fg;   tmp.v.z=  G;tmp.n=n++;     vertex_list.push_back(tmp);//4
   tmp.v.x=  0;   tmp.v.y=-Fg;   tmp.v.z=  G;tmp.n=n++;     vertex_list.push_back(tmp);//5
   tmp.v.x=  0;   tmp.v.y=-Fg;   tmp.v.z= -G;tmp.n=n++;     vertex_list.push_back(tmp);//6
   tmp.v.x=  0;   tmp.v.y= Fg;   tmp.v.z= -G;tmp.n=n++;     vertex_list.push_back(tmp);//7
   tmp.v.x=  G;   tmp.v.y=  0;   tmp.v.z= Fg;tmp.n=n++;     vertex_list.push_back(tmp);//8
   tmp.v.x=  G;   tmp.v.y=  0;   tmp.v.z=-Fg;tmp.n=n++;     vertex_list.push_back(tmp);//9
   tmp.v.x= -G;   tmp.v.y=  0;   tmp.v.z=-Fg;tmp.n=n++;     vertex_list.push_back(tmp);//10
   tmp.v.x= -G;   tmp.v.y=  0;   tmp.v.z= Fg;tmp.n=n++;     vertex_list.push_back(tmp);//11
   //each edge length is 2.0
   #define RI_Edge_Length 2.0
   std::vector<std::pair<int,int> > RI_edges;
   for(size_t fst=0,fst_size=vertex_list.size();fst<fst_size;++fst)
   {
      for(size_t snd=0,snd_size=vertex_list.size();snd<snd_size;++snd)
      {
         if(fst<snd)
         {
            const double dis = (vertex_list[fst].v-vertex_list[snd].v).norm();
            if(fabs(dis-RI_Edge_Length)<D_ERROR)
            {
               RI_edges.push_back(
                  std::pair<int,int> (vertex_list[fst].n, vertex_list[snd].n)
               );
            }
         }
      }
   }
   int t_n=0;
   for(size_t e1=0,e1_size=RI_edges.size();e1<e1_size;++e1)
   {
      for(size_t e2=0,e2_size=RI_edges.size();e2<e2_size;++e2)
      {
         const std::pair<int,int>& A = RI_edges[e1];
         const std::pair<int,int>& B = RI_edges[e2];
         if((A.first!=B.first)&&(A.second!=B.second))
         {
            if(A.first==B.second)
            {
               const double dis = (vertex_list[A.second].v-vertex_list[B.first].v).norm();
               if(fabs(dis-RI_Edge_Length)<D_ERROR)
               {
                  Triangle tmp;
                  tmp.n=t_n;
                  tmp.a=A.first;
                  tmp.b=A.second;
                  tmp.c=B.first;
                  triangle_list.push_back(tmp);
               }
            }
         }
      }
   }
}

inline int Existence(const std::vector<Coordinate>& list, const Coordinate& cor)
{
   for(size_t i=0,size=list.size();i<size;++i)
   {
      if((list[i].v-cor.v).norm()<D_ERROR)
      {
         return list[i].n;
      }
   }
   return cor.n;
}


void DevideT(std::vector<Coordinate>& vertex_list,std::vector<Triangle>& triangle_list)
{
   const std::vector<Coordinate> v_master=vertex_list;
   const std::vector<Triangle>   t_master=triangle_list;
   std::vector<Coordinate> v_tmp;
   std::vector<Triangle>   t_tmp;
   int num=vertex_list.size();

   for(size_t i=0,size=triangle_list.size();i<size;++i)
   {
      const Triangle& pk = triangle_list[i];
      const Coordinate& a = vertex_list[pk.a];
      const Coordinate& b = vertex_list[pk.b];
      const Coordinate& c = vertex_list[pk.c];
      
      //make vertex  
      Coordinate ab; ab.v = (a.v+b.v)/2.0; ab.n=num++;
      const int jud_ab = Existence(v_tmp,ab);
      if(jud_ab==ab.n)
      {
         v_tmp.push_back(ab);
      }
      else{num--;ab.n=jud_ab;}

      Coordinate bc; bc.v = (b.v+c.v)/2.0; bc.n=num++;
      const int jud_bc = Existence(v_tmp,bc);
      if(jud_bc==bc.n)
      {
         v_tmp.push_back(bc);
      }
      else{num--;bc.n=jud_bc;}

      Coordinate ca; ca.v = (c.v+a.v)/2.0; ca.n=num++;
      const int jud_ca = Existence(v_tmp,ca);
      if(jud_ca==ca.n)
      {
         v_tmp.push_back(ca);
      }
      else{num--;ca.n=jud_ca;}

      //devide
      Triangle t1,t2,t3,t4;
      const Triangle master=pk;
      t1.a=master.a; t1.b=ab.n;     t1.c=ca.n;
      t2.a=ab.n;     t2.b=master.b; t2.c=bc.n;
      t3.a=ca.n;     t3.b=bc.n;     t3.c=master.c;
      t4.a=ab.n;     t4.b=bc.n;     t4.c=ca.n;
      t_tmp.push_back(t1);
      t_tmp.push_back(t2);
      t_tmp.push_back(t3);
      t_tmp.push_back(t4);
   }

   vertex_list.clear();
   vertex_list.insert(vertex_list.end(),v_master.begin(),v_master.end());
   vertex_list.insert(vertex_list.end(),v_tmp.begin(),v_tmp.end());
   triangle_list.clear();
   triangle_list.insert(triangle_list.end(),t_tmp.begin(),t_tmp.end());
}

std::pair<std::vector<Coordinate>,std::vector<Triangle> >
Sphere(const double radius)
{  //Template of Sphere
   std::vector<Coordinate> vertex_list;
   std::vector<Triangle> triangle_list;

   RI(vertex_list,triangle_list);
   //Split 3 times
   //for(int i=0;i<3;++i)
   for(int i=0;i<Devide_N;++i)
   {
      DevideT(vertex_list,triangle_list);
   }
   //Projection
   for(size_t i=0,size=vertex_list.size();i<size;++i)
   {
      vertex_list[i].v/=vertex_list[i].v.norm();
      vertex_list[i].v*=radius;
   }

   return std::pair<std::vector<Coordinate>,std::vector<Triangle> > (vertex_list,triangle_list);
}

//************************ tetra bond ************************************************
//std::pair<std::vector<Coordinate>, std::vector<Rectangle<int>> >
//std::vector<Coordinate>
std::vector<Vector3D>
make_Tube
(const Vector3D& a, const Vector3D& b, double Diameter)
{  
   //construction of bottom disk
   const double Delta_Angle = 2.0*M_PI/(double)Devide_tube;
   //std::vector<Coordinate> disk(Devide_tube);
   std::vector<Vector3D> disk(Devide_tube);
   Vector3D axis_ab=a-b;
   
   //First, to make first point that is rotated around ab axis.
   const double r = Diameter/2.0;
   double theta=angle(axis_ab,Vector3D(0,1,0));
   if(isnan(theta)){ theta=0.0001;  axis_ab=Vector3D(0.0001, 0.0, 0.0);}
   const Vector3D normal_ab_axis = axis_ab/axis_ab.norm();          
   const Vector3D point_on_xy((r*cos(theta)), r*sin(theta), 0.0);  
   const Vector3D rota_axis = (normal_ab_axis*r%point_on_xy)/(normal_ab_axis*r%point_on_xy).norm(); 
   const Quaternion q0(rota_axis, M_PI/2.0);
   const Quaternion first_rota_p= q0*(normal_ab_axis*r)*bar(q0);
   //disk[0].v = Vector3D(first_rota_p.b, first_rota_p.c, first_rota_p.d);
   disk[0] = Vector3D(first_rota_p.b, first_rota_p.c, first_rota_p.d);

   //Second, to make another point that is rotated around ab axis.
   const Quaternion q(normal_ab_axis, Delta_Angle);
   for(int i=1; i<Devide_tube; ++i)
   {
      //const Quaternion result=q*disk[i-1].v*bar(q);
      const Quaternion result=q*disk[i-1]*bar(q);
      //disk[i].v=Vector3D(result.b, result.c, result.d);
      disk[i]=Vector3D(result.b, result.c, result.d);
   }

   auto top    = disk;
   auto bottom = disk;
   //int i = 0;
   //std::for_each(top.begin(),top.end(),      [&](Coordinate& h){h.v += a; h.n= i++;});
   std::for_each(top.begin(),top.end(),      [&](Vector3D& h){h += a;});
   //std::for_each(bottom.begin(),bottom.end(),[&](Coordinate& h){h.v += b; h.n= i++;});
   std::for_each(bottom.begin(),bottom.end(),[&](Vector3D& h){h += b;});

   //mesh
   //std::vector<Rectangle<int> > box(Devide_tube);
   //for(int i=0; i<Devide_tube; ++i)
   //{
   //   box[i%Devide_tube].a = top[i%Devide_tube].n;
   //   box[i%Devide_tube].b = top[(i+1)%Devide_tube].n;
   //   box[i%Devide_tube].c = bottom[i%Devide_tube].n;
   //   box[i%Devide_tube].d = bottom[(i+1)%Devide_tube].n;
   //}  
   top.insert(top.end(), bottom.begin(), bottom.end());
   //return std::pair<std::vector<Coordinate>, std::vector<Rectangle<int>> >(top, box); 
   //return std::vector<Coordinate>(top); 
   return std::vector<Vector3D>(top); 
}

std::pair<Vector3D,double> NearestOnPlane
(const Vector3D& point, const Vector3D& Pnormal, const Vector3D& on_plane_p)
{
   //normalized surface normal
   const Vector3D normal = Pnormal/Pnormal.norm();
   const double& a = normal.x;
   const double& b = normal.y;
   const double& c = normal.z;
   const double d = -(a*on_plane_p.x+b*on_plane_p.y+c*on_plane_p.z);

   const double t = -(a*point.x+b*point.y+c*point.z+d);
   const Vector3D vt(a*t,b*t,c*t);
   return std::pair<Vector3D, double> ((point+vt),vt.norm());
}


