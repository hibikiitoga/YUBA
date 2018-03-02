#ifndef TextReader_H
#define TextReader_H

#include "structure.hpp"

namespace TextReader
{
std::vector< std::string >
str_list_fast
(const std::string& file_name, const std::string& tag1, const std::string& tag2)
{
   std::vector< std::string > results;
   std::ifstream ifs(file_name);
   if(ifs.fail()){std::cout<<"file name fail"<<std::endl;}

   int  skip_line = 0    ; //the number of lines which data follows
   bool read_f    = false; //which data is target?
   bool fin_f     = false; //collection finished?

   std::string reading;    //a extracted line
   while(getline(ifs,reading)&&!fin_f)
   {
      if(skip_line==0)//tags or brank
      {
         std::vector< std::string > vsplit(4);
         boost::algorithm::split(vsplit,reading,boost::is_any_of(" ,\t"));
         if(vsplit.size()>=(size_t)4&&vsplit[0]=="#")
         {
               skip_line = boost::lexical_cast<int> (vsplit[3]);
            if(vsplit[1]==tag1&&vsplit[2]==tag2)
            {
               read_f = true;
            }else if(read_f){ fin_f=true;}
         }
      }
      else
      {
         if(read_f)
         {
            results.push_back(reading);
         }
         --skip_line;
      }
   }
   return results;
}


std::vector< std::string >
str_list_fast
(const std::string& file_name, const std::string& tag1, const int& step)
{
   return str_list_fast(file_name,tag1,boost::lexical_cast<std::string> (step));
}


std::vector< std::string >
str_list
(const std::string& file_name, const std::string& tag1, const std::string& tag2)
{
   //Extract partial data/info, tag1 is type of data

   std::vector< std::string > result;

   std::ifstream ifs(file_name);
   if(ifs.fail()){std::cout<<"file name fail"<<std::endl;}

   std::string temp;
   bool read_f = false;
   bool tags_f = false;
   while(getline(ifs,temp))
   {
      std::vector<std::string> vsplit;
      boost::algorithm::split(vsplit,temp,boost::is_any_of(" ,\t"));
      //read_f is "read or not" flag 
      if(vsplit.size()>=(size_t)3&&vsplit[0]=="#"&&vsplit[1]==tag1&&vsplit[2]==tag2)
      {
         read_f=true;
         tags_f=true;
      }
      else if(vsplit[0]=="#"){read_f=false;}
      else{tags_f=false;}
      
      if(read_f&&!tags_f&&!(vsplit.size()<(size_t)2))
      {
         result.push_back(temp); 
      }
   }

   return result;
}

std::vector< std::string >
str_list
(const std::string& file_name, const std::string& tag1, const int& step)
{
   return str_list(file_name,tag1,boost::lexical_cast<std::string> (step));
}


Coordinate cast_Coordinate(const std::string& req)
{
   Coordinate result;
   std::vector<std::string> vsplit;
   boost::algorithm::split(vsplit, req, boost::is_any_of(" ,\t"));
   result.n   = boost::lexical_cast<int>    (vsplit[0]);
   result.v.x = boost::lexical_cast<double> (vsplit[1]);
   result.v.y = boost::lexical_cast<double> (vsplit[2]);
   result.v.z = boost::lexical_cast<double> (vsplit[3]);
   return result;
}

Vector3D cast_Vector3D(const std::string& req)
{//text Coordinates type -> Vector3D 
   std::vector<std::string> vsplit;
   boost::algorithm::split(vsplit, req, boost::is_any_of(" ,\t"));
   return Vector3D
   (
   boost::lexical_cast<double> (vsplit[1]),
   boost::lexical_cast<double> (vsplit[2]),
   boost::lexical_cast<double> (vsplit[3])
   );
}

std::vector<Coordinate> cast_Coordinate(const std::vector<std::string>& req)
{
   std::vector<Coordinate> result;
   for(size_t i=0,size=req.size();i<size;++i)
   {
      result.push_back(cast_Coordinate(req[i]));
   }
   return result;
}

Triangle cast_Triangle(const std::string& req)
{
   Triangle result;
   std::vector<std::string> vsplit;
   boost::algorithm::split(vsplit, req, boost::is_any_of(" ,\t"));
   result.n  = boost::lexical_cast<int>   (vsplit[0]);
   result.a  = boost::lexical_cast<int>   (vsplit[1]);
   result.b  = boost::lexical_cast<int>   (vsplit[2]);
   result.c  = boost::lexical_cast<int>   (vsplit[3]);
   if(vsplit.size()==7)
   {
      result.AB = boost::lexical_cast<int>   (vsplit[4]);
      result.BC = boost::lexical_cast<int>   (vsplit[5]);
      result.CA = boost::lexical_cast<int>   (vsplit[6]);
   }
   return result; 
}

std::vector<Triangle> cast_Triangle(const std::vector<std::string>& req)
{
   std::vector<Triangle> result;
   for(size_t i=0,size=req.size();i<size;++i)
   {
      result.push_back(cast_Triangle(req[i]));
   }
   return result;
}


//************************ Tube **********************************************
Bond bond_pair 
(const std::string& req)
{
   Bond result;
   std::vector<std::string> vsplit;
   boost::algorithm::split(vsplit, req, boost::is_any_of(" ,\t"));
   result.v1 = boost::lexical_cast<int>   (vsplit[0]);
   result.v2 = boost::lexical_cast<int>  (vsplit[1]);
   if(scalar_f)
   {
      result.scalar=boost::lexical_cast<int> (vsplit[2]);
   }
   else
   {
      result.scalar=-1;
   }
   return result;
}

std::vector<Bond> cast_Bond
(const std::vector<std::string>& req)
{
   std::vector<Bond> bond_pair_list;

   for(size_t i=0, size=req.size(); i<size; ++i)
   {  
      Bond pair=bond_pair(req[i]);

      bool over=false;
      for(size_t j=0; j<bond_pair_list.size() && !over; ++j)
      {
         if((pair.v1==bond_pair_list[j].v2) 
            && (pair.v2==bond_pair_list[j].v1))
         {  over=true; }
      }
      if(!over)
      {
         bond_pair_list.push_back(pair);
      }
   }
   return bond_pair_list;
}
//cast scalar-----------------------------
int cast_sc_bead(const std::string& req)
{
   int result;

   std::vector<std::string> vsplit;
   boost::algorithm::split(vsplit, req, boost::is_any_of(" ,\t"));

   return result=boost::lexical_cast<int> (vsplit[4]);
}
std::vector<int> cast_sc_bead(const std::vector<std::string>& req)
{
   std::vector<int> result;
   for(size_t i=0, size=req.size(); i<size; ++i)
   {  
      result.push_back(cast_sc_bead(req[i]));
   }
   return result;
}

}//namespace TextReader end
#endif

