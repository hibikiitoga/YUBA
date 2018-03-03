#include <iostream>
#include <cmath>
#include <tuple>
#include <vector>
#include "Vector3Dedit.hpp"
#include "Quaternion.hpp"

//typedef struct
//{
//   int n;
//   int a;
//   int b;
//   int c;
//}Triangle;
typedef struct
{
   int n;
   int a;
   int b;
   int c;
   int d;
}Quadrangle;

std::tuple
<
   std::vector<Vector3D>, 
   std::vector<Triangle>, 
   std::vector<Quadrangle>
>
Regular_Sphere
(const int number_of_split)
{
   //make ling
   const std::vector<Vector3D> ling =
      [&number_of_split]()
      {
         const Vector3D Yaxis(0,1,0);
         const double delta_n = 2.0*M_PI/number_of_split;
         Vector3D front(1,0,0);
         Quaternion q(Yaxis,delta_n);
         std::vector<Vector3D> results;
         const auto q2v = [](const Quaternion& q_){return Vector3D(q_.b,q_.c,q_.d);};
         //for(int i=0;i<n;++i)
         for(int i=0;i<number_of_split;++i)
         {
            results.push_back(front);
            front=q2v(q*front*bar(q));
         }
         return results;
      }();
   const std::vector<double> middle_ys =
      [&number_of_split]()
      {
         std::vector<double> results;
         const double delta_n = 
            (number_of_split%2)?(2*M_PI/number_of_split):(2*M_PI/(number_of_split+1));
         double y;
         for(int i=1;i<number_of_split;++i)
         {
            y=std::abs(sin(delta_n*i));
            results.push_back(+y);
            results.push_back(-y);
         }
         results.push_back(0.0);
         std::sort(results.begin(),results.end(),[](const double& a, const double& b){return a>b;});
         return results;
      }();
   const std::vector<Vector3D> points = 
      [&]()
      {
         std::vector<Vector3D> results; 
         results.push_back(Vector3D(0,+1,0));
         for(int i=0,size=middle_ys.size();i<size;++i)
         {
            const double& y = middle_ys[i];
            const double factor = std::sqrt(1-y*y);
            for(int l=0,slze=ling.size();l<slze;++l)
            {
               results.push_back(Vector3D(factor*ling[l].x,y,factor*ling[l].z));
            }
         }
         results.push_back(Vector3D(0,-1,0));
         return results;
      }();
   const std::vector<Triangle> triangles =
      [&]()
      {
         std::vector<Triangle> results;
         const int top    = 0;
         const int bottm  = ((int)points.size()-1);
         Triangle tmp;
         int count=0;
         for(int i=0,size=ling.size();i<size;++i)
         {
            tmp.n=count++;
            tmp.a=top;
            tmp.b=(i)+1;
            tmp.c=((i+1)%size)+1;
            results.push_back(tmp);
         }
         for(int i=0,size=ling.size();i<size;++i)
         {
            tmp.n=count++;
            tmp.a=bottm;
            tmp.b=bottm-((i)+1);
            tmp.c=bottm-(((i+1)%size)+1);
            results.push_back(tmp);
         }

         return results;
      }();
   const std::vector<Quadrangle>  quadrangles =
      [&]()
      {
         std::vector<Quadrangle> results;
         Quadrangle tmp;
         int count=0;
         for(int i=0,size=middle_ys.size();i<(size-1);++i)
         {
            for(int l=0,slze=ling.size();l<slze;++l)
            {
               const int f=l;
               const int e=(l+1)%slze;
               tmp.n = count++;
               tmp.a = 1+f+i*slze;
               tmp.b = 1+e+i*slze;
               tmp.d = 1+f+(i+1)*slze;
               tmp.c = 1+e+(i+1)*slze;
               results.push_back(tmp);
            }
         }
         return results;
      }();

   return std::tuple<std::vector<Vector3D>,std::vector<Triangle>,std::vector<Quadrangle> > (points,triangles,quadrangles);
}

std::tuple
<
   std::vector<Vector3D>, 
   std::vector<Triangle>, 
   std::vector<Quadrangle>
>
Ellipsoid
(
   const double a,
   const double b,
   const int number_of_split
)
{
   std::tuple
   <
      std::vector<Vector3D>, 
      std::vector<Triangle>, 
      std::vector<Quadrangle>
   > result = Regular_Sphere(number_of_split*20);
  
   for(int i=0,size=std::get<0>(result).size();i<size;++i)
   {
      (std::get<0>(result)).at(i) *= a;
   }
   
   for(int i=0,size=std::get<0>(result).size();i<size;++i)
   {
      if(a>b)
      {
         (std::get<0>(result)).at(i).x *= b/a;
         (std::get<0>(result)).at(i).z *= b/a;
      }else{
         (std::get<0>(result)).at(i).x *= a/b;
         (std::get<0>(result)).at(i).z *= a/b;
      }
   }

   return result;
}
