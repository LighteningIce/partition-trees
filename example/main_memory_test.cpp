#include <iostream>
// #include <cstring>
#include <string>
using namespace std;
class String
{
private:
    char* _data;    // 存储字符串数据
    unsigned int _size;      // 保存字符串的长度
    unsigned int _capacity;  // 保存当前内存空间的容量
public:
    String();                           // 默认构造函数
    String(const char* str);            // 从字符数组创建字符串对象
    String(const String &other);        // 拷贝构造函数
    ~String();                          // 析构函数
    unsigned int length()const;                       // 返回字符串长度
    void print()const;                       // 打印字符串，末尾需要无换行
    void append(char ch);               // 添加字符
    String concat(const String &other)const; // 将自己与other合并后返回新的字符串
    String subString(int start, int length)const; // 从start开始截取长度为length的字符串
};

String::String()
{
 _capacity = 2 ;    _size = 0 ;
 _data = new char[_capacity] ;
 _data[_size] = '\0' ;
}

String::String(const char* str)
{
 _capacity = 2 ;    _size = 0 ;
 for(_size = 0 ; str[_size]!='\0' ; _size++);
 while(_capacity <= _size)
  _capacity <<= 1 ;
 _data = new char[_capacity] ;
 for(int i = 0 ; i<_size ; i++)
  _data[i] = str[i];
 _data[_size] = '\0' ;
}

String::String(const String &other)
{
 _capacity = other._capacity ;
 _size = other._size ;
 _data = new char[_capacity] ;
 for(int i = 0 ; i<_size ; i++)
  _data[i] = other._data[i] ;
 _data[_size] = '\0' ;
}

String::~String()
{
    delete []_data ;
}

unsigned int String::length()const
{
 return _size ;
}

void String::print()const
{
 for(int i = 0 ; i<_size ; i++)
  cout<<_data[i] ;
 return ;
}

void String::append(char ch)
{
 if(_capacity > _size+1)
 {
  _data[_size++] = ch ;
  _data[_size] = '\0' ;
 }
 else
 {
  while(_capacity <= _size+1)  _capacity <<= 1 ;
  char *newdata = new char[_capacity] ;
  for(int i = 0 ; i<_size ; i++)
  {
   newdata[i] = _data[i] ;
  }
  newdata[_size++] = ch ;
  newdata[_size] = '\0' ;
     delete []_data;
  _data = newdata ;
  newdata = NULL ;
 }
}

String String::concat(const String &other)const
{
    String *newstring = new String(_data) ;
    while(newstring->_capacity <= _size+other._size)   newstring->_capacity <<= 1 ;
    for(int i = 0; i<_size ; i++)
        newstring->_data[i] = _data[i] ;
    for(int i = 0; i<other._size ; i++)
        newstring->_data[newstring->_size++] = other._data[i] ;
    newstring->_data[newstring->_size] = '\0' ;
    return *newstring ;
}

String String::subString(int start, int length)const
{
    char *newdata = new char[length+1] ;
    int newsize = 0 ;
    for( ; newsize<length && (newsize+start)<_size ; newsize++)
        newdata[newsize] = _data[start+newsize] ;
    newdata[newsize] = '\0' ;
    String *newstring = new String(newdata) ;
    delete []newdata ;
    return *newstring ;
}


int main(){
    char str1[1001] = {0} ;
    char str2[1001] = {0} ;
    cin >> str1 >> str2 ;
    int M = strlen(str1);
    int N = strlen(str2) ;
    String s1(str1) ;
    String s2 ;
    for(int i = 0 ; i < N ; i++)
        s2.append(str2[i]) ;
    String s3 = s1.concat(s2) ;
    String s4 = s3.subString(M , N) ;
    s3.print() ;
    cout << endl ;
    s4.print() ;
    cout << endl ;
}