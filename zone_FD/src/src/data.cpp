#include "data.hpp"

/**
 * @brief 带参数的构造函数实现
 * @param n_ 网格点数量
 * @param nVar_ 变量数量
 */
Data::Data(int n_, int nVar_) {
  n = n_;
  nVar = nVar_;
  data.resize(nVar * n, 0.0);
}

/**
 * @brief 解初始化函数实现，用于设置解的初始条件
 * 
 * 根据不同的测试案例设置初始条件，例如Sod激波管问题
 * @param n_ 网格点数量
 * @param nvar_ 变量数量
 */
void Data::solInit(int n_, int nvar_) {
  n = n_;
  nVar = nvar_;
  data.resize(nVar * n, 0.0);

  for (int i = 0; i < n; i++) {
    real h = 2.0 / n;
    real xi = h / 2.0 + i * h - 1.0;
    // data[i]=-sin(M_PI*xi);//for burgers equation

    // for sod tube 1D

    real gamma = GAMMA;
    if (xi < 0) {
      (*this)(i, 0) = 1.0;
      (*this)(i, 1) = 0;
      (*this)(i, 2) = 1.0 / (gamma - 1) * 1;
    } else {
      (*this)(i, 0) = 0.125;
      (*this)(i, 1) = 0;
      (*this)(i, 2) = 1.0 / (gamma - 1) * 0.1;
    }
  }
}

/**
 * @brief 初始化函数实现
 * @param n_ 网格点数量
 * @param nvar_ 变量数量
 */
void Data::init(int n_, int nvar_) {
  n = n_;
  nVar = nvar_;
  data.resize(nVar * n, 0.0);
}

/**
 * @brief 重载括号操作符实现，用于访问特定网格点和变量的数据
 * @param i 网格点索引
 * @param ivar 变量索引
 * @return 对应位置数据的引用
 */
real &Data::operator()(int i, int ivar) { return data[i * nVar + ivar]; }

/**
 * @brief 重载方括号操作符实现，用于按线性索引访问数据
 * @param i 线性索引
 * @return 对应位置数据的引用
 */
real &Data::operator[](int i) {
  // return (this->data)[i];
  return data.at(i);
}

/**
 * @brief 计算指定变量的最大值实现
 * @param ivar 变量索引
 * @return 指定变量的最大值
 */
real Data::maxElement(int ivar) {
  if (data.empty())
    return 0;
  real res = data[0];
  for (int i = 1; i < n; i++) {
    res = (res < data[i]) ? data[i] : res;
  }
  return res;
}

/**
 * @brief 设置数据值实现
 * @param value 包含所有数据的向量
 */
void Data::setValue(std::vector<real> value) {
  assert(value.size() == data.size());
  std::copy(value.begin(), value.end(), data.begin());
}

// void Data::operator= (Data& dat)
// {
//     if(this->n==dat.n&&this->nVar==dat.nVar)
//     {
//         for(ind i = 0; i < n*nVar; i++)
//         {
//             (*this)[i]=dat[i];
//         }

//     }
//     else
//     {
//         std::cout<<"incorrect size at class Data=Data\n";
//     }
// }

/**
 * @brief 重载加等于操作符实现，将向量数据加到当前数据上
 * @param arr 要加上的数据向量
 */
void Data::operator+=(std::vector<real> arr) {
  assert(arr.size() == n * nVar);
  for (int i = 0; i < n * nVar; i++) {
    (*this)[i] += arr[i];
  }
}

/**
 * @brief 将所有数据设置为零的实现
 */
void Data::setZeros() { std::fill(data.begin(), data.end(), 0.0); }

/**
 * @brief 获取数据起始迭代器实现
 * @return 指向第一个数据元素的迭代器
 */
std::vector<real>::iterator Data::begin() { return data.begin(); }

/**
 * @brief 获取指定网格点数据的迭代器实现
 * @param i 网格点索引
 * @return 指向第i个网格点数据的迭代器
 */
std::vector<real>::iterator Data::random(int i) {
  return data.begin() + i * nVar;
}

/**
 * @brief 获取数据结束迭代器实现
 * @return 指向最后一个数据元素之后的迭代器
 */
std::vector<real>::iterator Data::end() { return data.end(); }

/**
 * @brief 获取数据总大小实现
 * @return 数据总大小（网格点数×变量数）
 */
int Data::size() { return data.size(); }

/**
 * @brief 获取变量数量实现
 * @return 变量数量
 */
int Data::getNVar() { return nVar; }

/**
 * @brief 获取网格点数量实现
 * @return 网格点数量
 */
int Data::getN() { return n; }

/**
 * @brief 提取特定变量的所有网格点数据实现
 * @param ivar 变量索引
 * @return 包含指定变量在所有网格点上数据的向量
 */
std::vector<real> Data::getIVar(int ivar) {
  std::vector<real> res;
  assert(ivar <= nVar);

  res.reserve(n);
  for (int i = 0; i < n; i++) {
    res.push_back((*this)(i, ivar));
  }
  return res;
}

/**
 * @brief 拷贝构造函数实现
 * @param origin_ 被拷贝的Data对象
 */
Data::Data(const Data &origin_) {
  n = origin_.n;
  nVar = origin_.nVar;
  data.resize(n * nVar);
  std::copy(origin_.data.begin(), origin_.data.end(), data.begin());
}

/**
 * @brief 设置变量名称实现
 * @param varName_ 变量名称的向量（右值引用）
 */
void Data::setvarName(std::vector<std::string> &&varName_) {
  assert(nVar == varName_.size());
  varName.resize(nVar);
  std::copy(varName_.begin(), varName_.end(), varName.begin());
}

/**
 * @brief 计算指定变量的L2范数实现
 * @param ivar 变量索引
 * @return 指定变量的L2范数
 */
real Data::getL2(int ivar) {
  if (ivar >= nVar) {
    std::cout << "Data error: ivar >= nVar";
    return 0;
  }
  real res = 0;
  for (int i = 0; i < n; i++) {
    res += pow((*this)(i, ivar), 2);
  }
  res = std::sqrt(res / n);
  return res;
}

/**
 * @brief 计算指定变量的L∞范数实现
 * @param ivar 变量索引
 * @return 指定变量的L∞范数
 */
real Data::getLinf(int ivar) {
  if (ivar >= nVar) {
    std::cout << "Data error: ivar >= nVar";
    return 0;
  }
  real res = std::abs((*this)(0, ivar));
  for (int i = 1; i < n; i++) {
    if (res < std::abs((*this)(i, ivar)))
      res = std::abs((*this)(i, ivar));
  }
  return res;
}