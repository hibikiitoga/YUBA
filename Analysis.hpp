#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP
#include <iostream>
#include <cmath>
#include <complex>
#include <tuple>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
typedef boost::multiprecision::cpp_int cpp_int;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50> > cdf;
#include "Vector3Dedit.hpp"
#include "Vector3Dcdf.hpp"
#include <boost/range/algorithm_ext/erase.hpp>
constexpr int maximum_n = 5;

namespace SD
{
   int SD_tau = INT_MAX;
   Getline sd_past_stream;
   Getline sd_curr_stream;
   std::pair<std::vector<unsigned int>,std::vector<unsigned int>> target_pair;//(s-f)^2;

   bool check();
   std::vector<double> calc_sd(const std::vector<Vector3D>& a,const std::vector<Vector3D>& b);

   bool check()
   {
      if(INT_MAX==SD_tau){std::cout<<"A"<<std::endl;return false;}
      if(!sd_past_stream.is_open()){std::cout<<"B"<<std::endl;return false;}
      if(!sd_curr_stream.is_open()){std::cout<<"C"<<std::endl;return false;}
      if(target_pair.first.empty()){std::cout<<"D"<<std::endl;return false;}
      return true;
   }

   void extract_steps
   (
      const std::set<unsigned int>& all,
      std::vector<unsigned int>& reqest
   )
   {
      std::vector<unsigned int> minus_tau_req = reqest;    
      std::for_each(minus_tau_req.begin(),minus_tau_req.end(),[](auto& r){r-=SD_tau;});
      const auto is_exist = [&](const unsigned int& t)->bool
      {
         for(auto i=all.begin();i!=all.end();++i)
         {
            if((*i)==t){return true;}
         }
         return false;
      };
      for(int i=0,size=reqest.size();i<size;++i)
      {
         if(is_exist(minus_tau_req[i]))
         {
            target_pair.first.push_back(minus_tau_req[i]);
            target_pair.second.push_back(reqest[i]);
         }
      }
   }

   std::vector<double> calc_sd(const std::vector<Vector3D>& a,const std::vector<Vector3D>& b)
   {
      if(a.size()!=b.size()){return std::vector<double>();}
      std::vector<double> result;
      result.reserve(a.size());
      for(int i=0,size=a.size();i<size;++i)
      {
         result.push_back((b[i]-a[i]).norm2());
      }
      return result;
   }
}

namespace SHF//define
{
   typedef Vector3Dcdf V3D;
   typedef Vector3D SphV;

   SphV inline rect2sphV(const Vector3D& rectV);
   template <typename T> inline constexpr int sign(const T& a);
   
   Vector3D inline Coordinate2Vector3D(const std::string& str);
   template<typename T>inline T cast(const std::string& str);
   template<>inline  Coordinate cast<Coordinate>(const std::string& str);
   template<>inline  Triangle cast<Triangle>(const std::string& str);

   cpp_int inline factorial(const int a);
   template<typename T> inline std::complex<T> exp_complex(const std::complex<T>& z);
   template<typename T> inline std::complex<T>
      multiplication
      (
         const std::complex<T>& cmp1, const std::complex<T>& cmp2
      );

   std::pair<std::vector<int>,std::vector<int> > inline
      ScanT
      (
         const int& pos,
         const std::vector<std::vector<int> >& vvNeighborPolygon,
         const std::vector<Triangle>& triangles
      );

   double inline
      sigma_i
      (
         const int& i,
         const std::vector<std::vector<int>>& vvNeighborPolygon,
         const std::vector<Triangle>& triangles,
         const std::vector<Vector3D>& vertice
      );

   std::complex<cdf> inline
      Y_fast
      (
         const int& n, const int& m,
         const cdf& getA, const cdf& getB, const cdf& getD,
         const cdf& theta, const cdf& phi
      );
   std::complex<cdf> inline 
      CCY_fast
      (
         const int& n, const int& m,
         const cdf& getA, const cdf& getB, const cdf& getD,
         const cdf& theta, const cdf& phi
      );


   std::vector<double> inline
      calc_us
      (const std::vector<Vector3D>& vertices,const std::vector<Triangle>& triangles);
}

namespace SHF//strage
{
   const std::vector<cdf> products = []()
   {
      //std::vector<cpp_int> result;
      std::vector<cdf> result;
      result.reserve(2*maximum_n+1);
      result[0]=1;
      for(int a=1,saze=2*maximum_n+1;a<saze;++a)
      {
         //result[a]=a*result[a-1];
         //result[a]=factorial(a);
         result[a]=cdf(factorial(a));
      }
      return result;
   }();
}

namespace SHF//function
{
   template <typename T>
   inline constexpr int sign(const T& a)
   {
      return (a>=0.0)?(+1.0):(-1.0);
   }

   SphV rect2sphV(const Vector3D& rectV)
   {
      const double& x = rectV.x;
      const double& y = rectV.y;
      const double& z = rectV.z;
      //const cdf r     = boost::multiprecision::sqrt(sqr(x)+sqr(y)+sqr(z));
      //const cdf theta = boost::multiprecision::acos(z/r);
      //const cdf phi   = sign(y)*boost::multiprecision::acos(x/(boost::multiprecision::sqrt(sqr(x)+sqr(y))));
      const double r     = std::sqrt(sqr(x)+sqr(y)+sqr(z));
      const double theta = std::acos(z/r); 
      const double phi   = sign(y)*std::acos(x/std::sqrt(sqr(x)+sqr(y)));
      return SphV(r,theta,phi);
   }

   Vector3D Coordinate2Vector3D(const std::string& str)
   {
      Vector3D result;
      std::vector<std::string> vsplit;
      boost::algorithm::split(vsplit, str, boost::is_any_of(" ,\t"));
      result.x = boost::lexical_cast<double> (vsplit[1]);
      result.y = boost::lexical_cast<double> (vsplit[2]);
      result.z = boost::lexical_cast<double> (vsplit[3]);
      return result;
   }

   template<>
   Coordinate cast<Coordinate>(const std::string& str)
   {
      Coordinate result;
      std::vector<std::string> vsplit;
      boost::algorithm::split(vsplit, str, boost::is_any_of(" ,\t"));
      result.n   = boost::lexical_cast<int>    (vsplit[0]);
      result.v.x = boost::lexical_cast<double> (vsplit[1]);
      result.v.y = boost::lexical_cast<double> (vsplit[2]);
      result.v.z = boost::lexical_cast<double> (vsplit[3]);
      return result;
   }
   
   template<>
   Triangle cast<Triangle>(const std::string& str)
   {
      Triangle result;
      std::vector<std::string> vsplit;
      boost::algorithm::split(vsplit, str, boost::is_any_of(" ,\t"));
      result.n  = boost::lexical_cast<int>   (vsplit[0]);
      result.a  = boost::lexical_cast<int>   (vsplit[1]);
      result.b  = boost::lexical_cast<int>   (vsplit[2]);
      result.c  = boost::lexical_cast<int>   (vsplit[3]);
      result.AB = boost::lexical_cast<int>   (vsplit[4]);
      result.BC = boost::lexical_cast<int>   (vsplit[5]);
      result.CA = boost::lexical_cast<int>   (vsplit[6]);
      return result;
   }

   cpp_int factorial(const int a)
   {
      cpp_int x = 1;
      for (int i = 1; i <= a; ++i)
      {
         x *= i;
      }
      return x;
   }

   template<typename T>
   std::complex<T>
   exp_complex(const std::complex<T>& z)
   {
      return std::complex<T>
         (
            boost::multiprecision::exp(z.real())*(boost::multiprecision::cos(z.imag())),
            boost::multiprecision::exp(z.real())*(boost::multiprecision::sin(z.imag()))
         );
   }

   template<typename T>
   std::complex<T>
   multiplication(const std::complex<T>& cmp1, const std::complex<T>& cmp2)
   {
      return std::complex<T>
         (
            (cmp1.real()*cmp2.real()-cmp1.imag()*cmp2.imag()), 
            (cmp1.real()*cmp2.imag()+cmp2.real()*cmp1.imag())
         );
   }
   
   std::pair<std::vector<int>,std::vector<int> >
   ScanT
   (
      const int& pos,
      const std::vector<std::vector<int> >& vvNeighborPolygon,
      const std::vector<Triangle>& triangles
   )
   {
      //Choise of Neighbors
      const std::vector<int>& Triangle_set = vvNeighborPolygon[pos];
   
      //construct of results array
      std::vector<int> ordered_triangle_set;
      std::vector<int> ordered_bond_set;
      ordered_triangle_set.reserve(Triangle_set.size()+2);
      ordered_bond_set.reserve(Triangle_set.size()+1);
   
      //first setting
      const Triangle& prev = triangles[Triangle_set[0]];
      Triangle read  ;
      int      next ;
      int      pair ;
   
      //pair vertex selection
      if(pos==prev.a)
      {//pair is b/c
         ordered_bond_set.push_back(prev.b);//-1 edge
         ordered_bond_set.push_back(prev.c);//front edge
         pair = prev.c;
         next = prev.CA;
      }
      else if(pos==prev.b)
      {//pair is c/a
         ordered_bond_set.push_back(prev.c);
         ordered_bond_set.push_back(prev.a);
         pair = prev.a;
         next = prev.AB;
      }
      else
      {//pair is a/b
         ordered_bond_set.push_back(prev.a);
         ordered_bond_set.push_back(prev.b);
         pair = prev.b;
         next = prev.BC;
      }
      ordered_triangle_set.push_back(prev.n);
   
      do{//Scan
         read = triangles[next];
         ordered_triangle_set.push_back(next);
         if(pos==read.a)
         {//new pair is b/c
            if(pair==read.b)
            {  pair=read.c; next=read.CA;   }
            else
            {  pair=read.b; next=read.AB;   }
         }
         else if(pos==read.b)
         {//new pair is c/a
            if(pair==read.c)
            {  pair=read.a; next=read.AB;   }
            else
            {  pair=read.c; next=read.BC;   }
         }
         else
         {//new pair is a/b
            if(pair==read.a)
            {  pair=read.b; next=read.BC;   }
            else
            {  pair=read.a; next=read.CA;   }
         }
         ordered_bond_set.push_back(pair);
      }while(next!=Triangle_set[0]);
      ordered_triangle_set.push_back(ordered_triangle_set[0]);
      ordered_bond_set.push_back(ordered_bond_set[1]);
   
      return std::pair<std::vector<int>,std::vector<int> > (ordered_triangle_set,ordered_bond_set);
   }
   
   double sigma_i
   (
      const int& i,
      const std::vector<std::vector<int>>& vvNeighborPolygon,
      const std::vector<Triangle>& triangles,
      const std::vector<Vector3D>& vertice
   )
   {
      const std::pair<std::vector<int>,std::vector<int> >& ScanT_result = ScanT(i,vvNeighborPolygon,triangles);
      const std::vector<int>& triangle_set = ScanT_result.first;
      const std::vector<int>&     edge_set = ScanT_result.second;
   
      double result = 0.0;
   
      for(size_t t=0,stze=triangle_set.size()-1;t<stze;++t)
      {
         const int& pos_i     = i;
         const int& pos_k     = edge_set[t];
         const int& pos_j     = edge_set[t+1];
         const int& pos_kdash = edge_set[t+2];
   
         //Substitute Coordinate
         const Vector3D vertex_k       = vertice[pos_k];
         const Vector3D vertex_i       = vertice[pos_i];
         const Vector3D vertex_j       = vertice[pos_j];
         const Vector3D vertex_kdash   = vertice[pos_kdash];
         const double&  l_ij           = (vertex_i-vertex_j).norm();
         const double&  theta1         = angle((vertex_i-vertex_k),(vertex_j-vertex_k));
         const double&  theta2         = angle((vertex_i-vertex_kdash),(vertex_j-vertex_kdash));
         const double&  sigma_ij       = l_ij*((1.0/std::tan(theta1))+(1.0/std::tan(theta2)))*0.5;
                        result        += (sigma_ij*l_ij)*0.25;
      }
      return result;
   }
   
   std::complex<cdf>
   Y_fast
   (
      const int& n,
      const int& m,
      const cdf& getA,
      const cdf& getB,
      const cdf& getD,
      const cdf& theta,
      const cdf& phi
   )
   {
      const cdf t = boost::multiprecision::cos(theta);
      cdf sum=0.0;
      for(int j=0;j<=(n-std::abs(m))/2;++j)//floor function, j<=max{n = (k-m)/2,n in Z}
      {
         const cdf s1(std::pow(-1,j));
         const cdf s2(products[(2*n-2*j)]);
         const cdf m1(products[j]);
         const cdf m2(products[(n-j)]);
         const cdf m3(products[(n-2*j-std::abs(m))]);
         const cdf powed  = boost::multiprecision::pow(t,n-2*j-std::abs(m)); 
         sum += ((s1*s2)/(m1*m2*m3))*powed;
      }
   
      const cdf Pk = getD*boost::multiprecision::pow((1.0-sqr(t)),std::abs(m)/2.0)*sum;
      const std::complex<cdf> c_i(0.0,1.0);
      const std::complex<cdf> c_m_phi(m*phi,0.0);
      const std::complex<cdf> exp_part = exp_complex(multiplication(c_i,c_m_phi));
      const cdf factor = getA*getB*Pk;
      const std::complex<cdf> result(factor*exp_part.real(),factor*exp_part.imag());
      return result;
   }

   std::complex<cdf>
   CCY_fast
   (
      const int& n,
      const int& m,
      const cdf& getA,
      const cdf& getB,
      const cdf& getD,
      const cdf& theta,
      const cdf& phi
   )
   {
      return std::conj(Y_fast(n,m,getA,getB,getD,theta,phi));
   }

   namespace source
   {//why cant stock! fuckking memory size!
      inline cdf get_A(const int& m);
      cdf get_A(const int& m)
      {//for extracting a part of (-1)^(m+|m|)/2
         return cdf(pow(-1.0,(double)(+m+std::abs(m))/2.0));
      }
   
      inline cdf get_B(const int& n,const int& m);
      cdf get_B(const int& n,const int& m)
      {
         const cdf factor = ((2.0*n+1.0)/(4.0*M_PI));
         const cdf sister = (factorial(n-std::abs(m))).convert_to<cdf>();
         const cdf mother = (factorial(n+std::abs(m))).convert_to<cdf>();
         return boost::multiprecision::sqrt((factor*sister)/mother);
      }
   
      inline cdf get_D(const int& n);
      cdf get_D(const int& n)
      {
         //return D[n];
         const cdf two(2.0);
         const cdf one(1.0);
         const cdf pow_cdf(boost::multiprecision::pow(two,n));
         return one/pow_cdf;
      }
   
      //const std::vector<cpp_int> products = []()
      const std::vector<cdf> products = []()
      {
         //std::vector<cpp_int> result;
         std::vector<cdf> result;
         result.reserve(2*maximum_n+1);
         result[0]=1;
         for(int a=1,saze=2*maximum_n+1;a<saze;++a)
         {
            //result[a]=a*result[a-1];
            //result[a]=factorial(a);
            result[a]=cdf(factorial(a));
         }
         return result;
      }();
   }//namespace source end
}

#endif
