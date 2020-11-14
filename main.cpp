#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <string>
#include <math.h>
#include <algorithm>
#include<ctime>
#include <stdio.h>
#include <sys/times.h>
#include <unistd.h>
#include <stdlib.h>
#include <map>
#include <set>
#include <malloc.h>

#define _GLIBCXX_USE_CXX11_ABI 0
#define n_real 4938921

using namespace std;

bool compare(pair<int,string> p1,pair<int,string> p2);
int main()
{
    ifstream f;
    string tmp;
    static char whole_gene[n_real];
    f.open("NC_008253.fna",ios::in);
    //f.open("testtext.fna",ios::in);

    if(f.fail()){
        cout << "can not open the file." << endl;
        exit(1);
    }
    //去掉第一行
    if(f.good())
        getline(f,tmp);
    int z=0;
    do{
        if(f.good()){
            getline(f,tmp);
            for(string::iterator it=tmp.begin();it<tmp.end();it++){
                whole_gene[z]=(*it);
                z++;
            }
        }
    }while(!f.eof());
    f.close();
    whole_gene[n_real-1]='$';
    //cout << whole_gene.size() << endl;

    int n=n_real;//加上最后的$
    int k=(int)20*log2(n);//组数
    int l=n/k;//每组的字符数
//    int l=6;
//    int k=4;

    int last_num = n-k*l;

    //给最后一组n-l(k-1) : n排序并计算CSA
    vector<pair<int,string> > T;
    int Ti[l+1];
    string s("");
    int sa = last_num-1;
    for(int i=n-1;i>=n-last_num;i--){
    	s.insert(s.begin(),whole_gene[i]);
    	cout << s << endl;
        T.push_back(make_pair(sa,s));
        sa--;
    }
    sort(T.begin(),T.end(),compare);
    int tl = T.size();

    cout << "l = " << l << endl;
    cout << "k = " << k << endl;
    cout << "num = " << last_num << endl;
    cout << "n = " << n << endl;


//    cout << "SA" << endl;
//    for(int i=0;i<tl;i++){
//        cout << T[i].first << " ";
//    }
//    cout << endl;
    //求SA-1
    int RSA[tl];//SA-1
    for(int i=0;i<tl;i++){
        RSA[T[i].first] = i;
    }
//    cout << "SA-1" << endl;
//    for(int i=0;i<tl;i++){
//        cout << RSA[i] << " ";
//    }
//    cout << endl;
    //求CSA值
    static int CSA[n_real];
    static int CSA_new[n_real];
    static int order_num[n_real] = {0};
    for(int i=0;i<tl;i++){
        int a = T[i].first;
        if(a<tl-1)
            CSA[i] = RSA[a+1];
        else
            CSA[i] = RSA[0];
    }
//    cout << "CSA" << endl;
//    for(int i=0;i<tl;i++){
//        cout << CSA[i] << " ";
//    }
//    cout << endl;

    //求分界点
    int depart_pos[4] = {-1};//用-1来标记
    depart_pos[0] = 0;
    for(vector<pair<int,string>>::iterator iter=T.begin();iter!=T.end();iter++){
        if((*iter).second.at(0)=='A')
            depart_pos[1] = iter - T.begin();
        if((*iter).second.at(0)=='C')
            depart_pos[2] = iter - T.begin();
        if((*iter).second.at(0)=='G')
            depart_pos[3] = iter - T.begin();
    }
    vector<pair<int,string>>().swap(T);


    int clktck = 0;
    clktck = sysconf(_SC_CLK_TCK);
    struct tms  tmsstart, tmsend;
    clock_t     Start, End;
    //merge step
    int order[l+1];//最后一个存的是T‘的order
    order[l] = RSA[0];
    //cout << order[l] << endl;
    //vector<pair<int,string> > Ti;
    for(int i=0;i<k;i++){//对每一组，从后往前编号

//        cout << "depart_pos" << endl;
//        for(int i=0;i<4;i++)
//            cout << depart_pos[i] << " ";
//        cout << endl;

        int basic = n - i*l-last_num-1;//每一组最后一个字符在whole_gene中的位置
        int basic_start = basic - l + 1;
        cout << "handling group "<<i+1<<endl;
        int T_size = i*l + last_num;//用T表示已经排好序的部分
        int depart_pos_increase[4]={0};

        Start = times(&tmsstart);
        //对该组中的每一个后缀，从后往前计算order值
         for(int j=l;j>0;j--){
            char c = whole_gene[basic-l+j];
            int st=order[j];
            int left,r,mid;
            switch(c){
                case 'A':
                    depart_pos_increase[1]++;
                    //depart_pos_increase[2]++;
                    //depart_pos_increase[3]++;
                    if(depart_pos[1]==-1 || CSA[depart_pos[0]+1]>order[j]){//B=空集
                        order[j-1]= depart_pos[0];
                        break;
                    }
//                    for(int k=depart_pos[0]+1;k<=depart_pos[1];k++){
//                        if(CSA[k]<=order[j]){
//                            order[j-1] = k;
//                        }
//                    }

                    left=depart_pos[0]+1;
                    r=depart_pos[1];
                    do{
                        mid = (left+r)/2;
                        if(CSA[mid]<=st)left=mid;
                        else if(CSA[mid]>st)r=mid;
                    }while(left!=r && left!=r-1);
                    if(st>=CSA[r])order[j-1] = r;
                    else order[j-1] = left;

                    break;
                case 'C':
                    depart_pos_increase[2]++;
                    //depart_pos_increase[3]++;
                    if(depart_pos[2]==-1 || CSA[depart_pos[1]+1]>order[j]){//B=空集
                        order[j-1]= depart_pos[1];
                        break;
                    }
//                    for(int k=depart_pos[1]+1;k<=depart_pos[2];k++){
//                        if(CSA[k]<=order[j]){
//                            order[j-1] = k;
//                        }
//
                    left=depart_pos[1]+1;
                    r=depart_pos[2];
                    do{
                        mid = (left+r)/2;
                        if(CSA[mid]<=st)left=mid;
                        else if(CSA[mid]>st)r=mid;
                    }while(left!=r && left!=r-1);
                    if(st>=CSA[r])order[j-1] = r;
                    else order[j-1] = left;
                    break;
                case 'G':
                    depart_pos_increase[3]++;
                    if(depart_pos[3]==-1 || CSA[depart_pos[2]+1]>order[j]){//B=空集
                        order[j-1]= depart_pos[2];
                        break;
                    }
//                    for(int k=depart_pos[2]+1;k<=depart_pos[3];k++){
//                        if(CSA[k]<=order[j]){
//                            order[j-1] = k;
//                        }
//                    }
                    left=depart_pos[2]+1;
                    r=depart_pos[3];
                    do{
                        mid = (left+r)/2;
                        if(CSA[mid]<=st)left=mid;
                        else if(CSA[mid]>st)r=mid;
                    }while(left!=r && left!=r-1);
                    if(st>=CSA[r])order[j-1] = r;
                    else order[j-1] = left;
                    break;
                case 'T':
                    depart_pos_increase[0]++;
                    if(depart_pos[3]+1==T_size){//之前没有以T开头的后缀
                        order[j-1]= depart_pos[3];
                        break;
                    }
                    if(CSA[depart_pos[3]+1]>order[j]){//B=空集
                        order[j-1]= depart_pos[3];
                        break;
                    }
//                    for(int k=depart_pos[3]+1;k<=i*l+last_num-1;k++){
//                        if(CSA[k]<=order[j]){
//                            order[j-1] = k;
//                        }
//                    }
                    left=depart_pos[3]+1;
                    r=i*l+last_num-1;
                    do{
                        mid = (left+r)/2;
                        if(CSA[mid]<=st)left=mid;
                        else if(CSA[mid]>st)r=mid;
                    }while(left!=r && left!=r-1);
                    if(st>=CSA[r])order[j-1] = r;
                    else order[j-1] = left;
                    break;
            }
        }

//        cout << "depart_pos_increase" << endl;
//        for(int i=0;i<4;i++)
//            cout << depart_pos_increase[i] << " ";
//        cout << endl;
//
//        cout << "order" << endl;
//        for(int z=0;z<l;z++){
//            cout << order[z] << "  ";
//        }
//        cout << endl;
        End = times(&tmsend);
        cout << "order time: " << (End - Start)/(double) clktck << endl;


        //建立f(j)原来T中现在的值
        //遍历order，计算order小于各个值的数量
        map<int,set<int>> order_map;
//        for(int x=0;z<l;z++){
//                set<int> st;
//            order_map.insert(x,st);
//        }
        Start = times(&tmsstart);
        for(int x=0;x<n_real;x++){
            order_num[x]=0;
        }
        for(int x=0;x<l;x++){
            order_num[order[x]]++;
            order_map[order[x]].insert(x);
        }
        int incr = order_num[0];
        order_num[0] = 0;
        for(int y=0;y<T_size-1;y++){
            int tmp = order_num[y+1];
            order_num[y+1]=order_num[y] + incr;
            incr = tmp;
        }
//        cout << "order_num" << endl;
//        int num_test=0;
//        for(int z=0;z<T_size;z++){
//            cout << order_num[z] << endl;
//            num_test += order_num[z];
//        }
//        cout << "num_test = " << num_test << endl;


//        cout << "f" << endl;
        int f[T_size];
        for(int k=0;k<T_size;k++){
            f[k] = k + order_num[k];
            //cout << f[k] << " ";
        }
//        cout << endl;
        End = times(&tmsend);
        cout << "f time: " << (End - Start)/(double) clktck << endl;

        //建立g(j)新加的后缀的值
        Start = times(&tmsstart);
        int g[l+1];
//        cout << "g" << endl;
//        int sa = l;//该组中最后一个字符的位置0-l
//        pair<int,string> Ti[l+1];
//        int d = 0;
//        for(int i=basic;i>basic-l;i--){
//            s.insert(s.begin(),whole_gene[i]);
////            cout << s << endl;
//            Ti[d] = make_pair(sa,s);
//            sa--;
//            d++;
//        }
//        sort(Ti,Ti+l+1,compare);//排的是原来本组新加的l个后缀
//        for(int k=1;k<=l;k++){
//            g[Ti[k-1].first] = order[Ti[k-1].first-1] + k;
//        }


        map<int,set<int>>::iterator iter;
        for(iter=order_map.begin();iter!=order_map.end();iter++){//对每个order值对应的后缀集合排序
            pair<int,string> r_set[1000];
            set<int> s = iter->second;
            set<int>::iterator iter_set;
            int m=0;
            for(iter_set=s.begin();iter_set!=s.end();iter_set++){
                char *s = (char*)whole_gene+basic_start+(*iter_set);
                if(s==NULL) cout << "no" <<endl;
                string str(s);
                //cout << str << endl;
                r_set[m]=make_pair(*iter_set,str);//x为该后缀在该分组中起始位置
                m++;
            }
            sort(r_set,r_set+m,compare);
            for(int k=0;k<m;k++){
                int r = order_num[iter->first]+k;
                Ti[r] = r_set[k].first+1;//Ti是该组中排序后的SA值
                g[r_set[k].first+1] = order[r_set[k].first] + r + 1;//1是加上自己
            }
        }


//        for(int x=1;x<=l;x++){
//            cout << g[x] << " ";
//        }
//        cout << endl;
        End = times(&tmsend);
        cout << "g time: " << (End - Start)/(double) clktck << endl;


        //计算新的CSA的值
        Start = times(&tmsstart);
        int T_size_new = T_size+l;
        //vector<int> CSA_new(T_size_new);
        int jf=1;int jg=1;
        for(int t=0;t<T_size_new;t++){
            if(t==0)
                CSA_new[t] = g[1];
            else if(t==g[l]){
                CSA_new[t] = f[CSA[0]];jg++;
            }
            else if(t==f[jf]){
                CSA_new[t] = f[CSA[jf]]; jf++;
            }
            else{
                CSA_new[t] = g[Ti[jg-1]+1];
                //cout << Ti[jg-1].first << " " << jg <<endl;
                jg++;
            }
        }
        End = times(&tmsend);
        cout << "update CSA time: " << (End - Start)/(double) clktck << endl;
        //CSA = CSA_new;
        for(int x=0;x<T_size_new;x++){
            CSA[x]=CSA_new[x];
        }

//        cout << "Ti" << endl;
//        for(int i=0;i<l;i++)
//            cout << Ti[i].first << " " << Ti[i].second;
//        cout << endl;
//
//        cout << "CSA" << endl;
//        for(int i=0;i<T_size_new;i++)
//            cout << CSA[i] << " ";
//        cout << endl;

        order[l] = g[1];

        //更新depart_pos
        for(int x=1;x<4;x++){
            if(depart_pos[x]==-1 && depart_pos_increase[x]!=0){
                depart_pos[x] = depart_pos[x-1];
//                for(int y=1;y<=x;y++){
//                    depart_pos[x] += depart_pos_increase[x];
//                }
                depart_pos[x] += depart_pos_increase[x];
            }
            else if(depart_pos[x]!=-1){
                for(int y=1;y<=x;y++){
                    depart_pos[x] += depart_pos_increase[y];
                }
            }
        }
        for(iter=order_map.begin();iter!=order_map.end();iter++){
                set<int>().swap(iter->second);
        }
        map<int,set<int>>().swap(order_map);
        //order_map.clear();
        malloc_trim(0);




    }

    ofstream outfile;
    outfile.open("result.txt");
    for(int i=0;i<n;i++)
        outfile << CSA[i] << " ";
    outfile.close();
    return 0;
}

bool compare(pair<int,string> p1, pair<int,string> p2){
    string s1 = p1.second;
    string s2 = p2.second;
	return s1+s2 < s2+s1;
}

