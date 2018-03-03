#ifndef Vector3D_CDF_H
#define Vector3D_CDF_H
#include<iostream>
#include<cmath>

typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50> > cdf;
cdf sqr(const cdf& x)
{
   return x*x;
}

//ベクトルの定義と各種計算
class Vector3Dcdf{
   public:
      cdf x;
      cdf y;
      cdf z;
         
      //コンストラクタ
      Vector3Dcdf();
      Vector3Dcdf(cdf x_,cdf y_,cdf z_);
         
      //代入演算子
      Vector3Dcdf& operator=(const Vector3Dcdf& v);
      
      //単項演算子
      Vector3Dcdf& operator+=(const Vector3Dcdf& v);
      Vector3Dcdf& operator-=(const Vector3Dcdf& v);
      Vector3Dcdf& operator*=(cdf k);
      Vector3Dcdf& operator/=(cdf k);
      Vector3Dcdf operator+()const;
      Vector3Dcdf operator-()const;

      //添字演算子
      double& operator[](int i);
      
      //比較演算子
      bool operator==(const Vector3Dcdf& v) const;
      bool operator!=(const Vector3Dcdf& v) const;

      //ベクトルの長さ
      cdf norm() const;
      cdf norm2()const;
      
      //正規化
      void normalize();
      //Vector3D normalize();
};

//ベクトル演算
Vector3Dcdf operator+(const Vector3Dcdf& u, const Vector3Dcdf& v);
Vector3Dcdf operator-(const Vector3Dcdf& u, const Vector3Dcdf& v);
Vector3Dcdf operator*(const cdf k, const Vector3Dcdf& v);
Vector3Dcdf operator*(const Vector3Dcdf& v, cdf  k);
Vector3Dcdf operator/(const Vector3Dcdf& v, cdf  k);
cdf operator*(const Vector3Dcdf& u, const Vector3Dcdf& v);
Vector3Dcdf operator%(const Vector3Dcdf& u, const Vector3Dcdf& v);
cdf angle(const Vector3Dcdf& u, const Vector3Dcdf& v);

//出力
inline std::ostream& operator<<(std::ostream& s, const Vector3Dcdf& v);

/********メンバ関数の実装**********************************************/
//コンストラクタ
inline Vector3Dcdf::Vector3Dcdf(){x=y=z=0.0;}
inline Vector3Dcdf::Vector3Dcdf(cdf x_, cdf y_, cdf z_){
   this->x=x_;  this->y=y_;  this->z=z_;
}
//代入演算子
inline Vector3Dcdf& Vector3Dcdf::operator=(const Vector3Dcdf& v){
   this->x=v.x;   this->y=v.y;   this->z=v.z;
   return *this;
}
//単項演算子
inline Vector3Dcdf& Vector3Dcdf::operator+=(const Vector3Dcdf& v){
   this->x += v.x;   this->y += v.y;   this->z += v.z;
   return *this;
}
inline Vector3Dcdf& Vector3Dcdf::operator-=(const Vector3Dcdf& v){
   this->x -= v.x;   this->y -= v.y;   this->z -= v.z;
   return *this;
}
inline Vector3Dcdf& Vector3Dcdf::operator*=(cdf k){
   this->x *= k;  this->y *= k;  this->z *= k;
   return *this;  
}
inline Vector3Dcdf& Vector3Dcdf::operator/=(cdf k){
   this->x /= k;  this->y /= k;  this->z /= k;
   return *this;  
}
inline Vector3Dcdf Vector3Dcdf::operator+()const{//+Vector3D
   return *this;
}
inline Vector3Dcdf Vector3Dcdf::operator-()const{//-Vector3D
   return Vector3Dcdf(-x,-y,-z);
}
////添字演算子
//inline double Vector3D::operator[](const int& i)const{
// if(i==0){
//          return x;
// }else if(i==1){
//          return y;
// }else if(i==2){
//          return z;
// }else{
//          exit(1);
// }
//}
//比較演算子
inline bool Vector3Dcdf::operator==(const Vector3Dcdf& v)const{
   return (x==v.x)&&(y==v.y)&&(z==v.z);
}
inline bool Vector3Dcdf::operator!=(const Vector3Dcdf& v)const{
   return !(*this==v);
}
//ベクトルの長さ
inline cdf Vector3Dcdf::norm()const{
   return boost::multiprecision::sqrt(sqr(x)+sqr(y)+sqr(z));
}
inline cdf Vector3Dcdf::norm2()const{
   return (sqr(x)+sqr(y)+sqr(z));
}
//正規化
inline void Vector3Dcdf::normalize(){
   *this /= norm();
}
/********グローバル関数の実装**********************************************/
//二項演算子の定義
//Vector3D+Vector3D
inline Vector3Dcdf operator+(const Vector3Dcdf& u, const Vector3Dcdf& v){
      Vector3Dcdf w;
      w.x=u.x+v.x;
      w.y=u.y+v.y;
      w.z=u.z+v.z;
      return w;
}
//Vector3D-Vector3D
inline Vector3Dcdf operator-(const Vector3Dcdf& u, const Vector3Dcdf& v){
      Vector3Dcdf w;
      w.x=u.x-v.x;
      w.y=u.y-v.y;
      w.z=u.z-v.z;
      return w;
}
//double*Vector3D
inline Vector3Dcdf operator*(cdf k, const Vector3Dcdf& v){
      return Vector3Dcdf(k*v.x, k*v.y, k*v.z);
}
//Vector3D*double
inline Vector3Dcdf operator*(const Vector3Dcdf& v, cdf k){
      return Vector3Dcdf(v.x*k, v.y*k, v.z*k);
}
//Vector3D/double
inline Vector3Dcdf operator/(const Vector3Dcdf& v,cdf k){
      return Vector3Dcdf(v.x/k, v.y/k, v.z/k);
}
//内積 Vector3D*Vector3D
inline cdf operator*(const Vector3Dcdf& u, const Vector3Dcdf& v){
      return u.x*v.x+u.y*v.y+u.z*v.z;
}
//外積 Vector3D%Vector3D
inline Vector3Dcdf operator%(const Vector3Dcdf& u, const Vector3Dcdf& v){
      Vector3Dcdf w;
      w.x=u.y*v.z-u.z*v.y;
      w.y=u.z*v.x-u.x*v.z;
      w.z=u.x*v.y-u.y*v.x;
      return w;
}
//画面への表示
inline std::ostream& operator<<(std::ostream& s, const Vector3Dcdf& v){
      return s<<v.x<<" "<<v.y<<" "<<v.z;
}
//２つのベクトルのなす角(rad)
inline cdf angle(const Vector3Dcdf& u, const Vector3Dcdf& v){
      cdf cos=u*v/(u.norm()*v.norm());
      return boost::multiprecision::acos(cos);
}
#endif
