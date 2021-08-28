
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#define random(a,b) (rand()%(b-a)+a)

using namespace std;
//g++ main_gen_data.cpp -o main_gen
//
int main(int argc, const char * argv[])
{   
    int num=atoi(argv[1]);
   unsigned int *array=new unsigned int[num+1];
   unsigned int *height=new unsigned int[num+1];
   for (int i = 0; i < num; i++)
   {
       array[i]=0;
       height[i]=0;
   }
   
   int max_degree=atoi(argv[2]);
   int max_height=atoi(argv[3]);
    srand((int)time(0));  // 产生随机种子  把0换成NULL也行
    cout <<1<<"   "<< 0<<"   "<<random(2, 8)<< "    "<<random(1,30)<<"    "<<random(1,50)<<endl;
    for (int i = 1; i < num; i++)
    {
        //id ,parent, nw, msw , ew
        int p;
        if(i==1)
        {
            p=1;
        }else
        {
            p=random(1, i) ;
        }
        while (1)
        {
            if(array[p]<=max_degree-1&&height[p]<=max_height-1){
                break;
            }else
            {
                p=random(1, i) ;
                
            }
            
        }
        
        array[p]++;
        height[i+1]=height[p]+1;
        cout <<i+1<<"   "<< p<<"   "<<random(2, 8)<< "    "<<random(1,200)<<"    "<<random(1,500)<<endl;
        sleep(0.5);
    }
    return 0;
}