<<<<<<< HEAD
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
    //å»æ‰ç¬¬ä¸€è¡Œ
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

    int n=n_real;//åŠ ä¸Šæœ€åçš„$
    int k=(int)20*log2(n);//ç»„æ•°
    int l=n/k;//æ¯ç»„çš„å­—ç¬¦æ•°
//    int l=6;
//    int k=4;

    int last_num = n-k*l;

    //ç»™æœ€åä¸€ç»„n-l(k-1) : næ’åºå¹¶è®¡ç®—CSA
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
    //æ±‚SA-1
    int RSA[tl];//SA-1
    for(int i=0;i<tl;i++){
        RSA[T[i].first] = i;
    }
//    cout << "SA-1" << endl;
//    for(int i=0;i<tl;i++){
//        cout << RSA[i] << " ";
//    }
//    cout << endl;
    //æ±‚CSAå€¼
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

    //æ±‚åˆ†ç•Œç‚¹
    int depart_pos[4] = {-1};//ç”¨-1æ¥æ ‡è®°
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
    int order[l+1];//æœ€åä¸€ä¸ªå­˜çš„æ˜¯Tâ€˜çš„order
    order[l] = RSA[0];
    //cout << order[l] << endl;
    //vector<pair<int,string> > Ti;
    for(int i=0;i<k;i++){//å¯¹æ¯ä¸€ç»„ï¼Œä»åå¾€å‰ç¼–å·

//        cout << "depart_pos" << endl;
//        for(int i=0;i<4;i++)
//            cout << depart_pos[i] << " ";
//        cout << endl;

        int basic = n - i*l-last_num-1;//æ¯ä¸€ç»„æœ€åä¸€ä¸ªå­—ç¬¦åœ¨whole_geneä¸­çš„ä½ç½®
        int basic_start = basic - l + 1;
        cout << "handling group "<<i+1<<endl;
        int T_size = i*l + last_num;//ç”¨Tè¡¨ç¤ºå·²ç»æ’å¥½åºçš„éƒ¨åˆ†
        int depart_pos_increase[4]={0};

        Start = times(&tmsstart);
        //å¯¹è¯¥ç»„ä¸­çš„æ¯ä¸€ä¸ªåç¼€ï¼Œä»åå¾€å‰è®¡ç®—orderå€¼
         for(int j=l;j>0;j--){
            char c = whole_gene[basic-l+j];
            int st=order[j];
            int left,r,mid;
            switch(c){
                case 'A':
                    depart_pos_increase[1]++;
                    //depart_pos_increase[2]++;
                    //depart_pos_increase[3]++;
                    if(depart_pos[1]==-1 || CSA[depart_pos[0]+1]>order[j]){//B=ç©ºé›†
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
                    if(depart_pos[2]==-1 || CSA[depart_pos[1]+1]>order[j]){//B=ç©ºé›†
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
                    if(depart_pos[3]==-1 || CSA[depart_pos[2]+1]>order[j]){//B=ç©ºé›†
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
                    if(depart_pos[3]+1==T_size){//ä¹‹å‰æ²¡æœ‰ä»¥Tå¼€å¤´çš„åç¼€
                        order[j-1]= depart_pos[3];
                        break;
                    }
                    if(CSA[depart_pos[3]+1]>order[j]){//B=ç©ºé›†
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


        //å»ºç«‹f(j)åŸæ¥Tä¸­ç°åœ¨çš„å€¼
        //éå†orderï¼Œè®¡ç®—orderå°äºå„ä¸ªå€¼çš„æ•°é‡
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

        //å»ºç«‹g(j)æ–°åŠ çš„åç¼€çš„å€¼
        Start = times(&tmsstart);
        int g[l+1];
//        cout << "g" << endl;
//        int sa = l;//è¯¥ç»„ä¸­æœ€åä¸€ä¸ªå­—ç¬¦çš„ä½ç½®0-l
//        pair<int,string> Ti[l+1];
//        int d = 0;
//        for(int i=basic;i>basic-l;i--){
//            s.insert(s.begin(),whole_gene[i]);
////            cout << s << endl;
//            Ti[d] = make_pair(sa,s);
//            sa--;
//            d++;
//        }
//        sort(Ti,Ti+l+1,compare);//æ’çš„æ˜¯åŸæ¥æœ¬ç»„æ–°åŠ çš„lä¸ªåç¼€
//        for(int k=1;k<=l;k++){
//            g[Ti[k-1].first] = order[Ti[k-1].first-1] + k;
//        }


        map<int,set<int>>::iterator iter;
        for(iter=order_map.begin();iter!=order_map.end();iter++){//å¯¹æ¯ä¸ªorderå€¼å¯¹åº”çš„åç¼€é›†åˆæ’åº
            pair<int,string> r_set[1000];
            set<int> s = iter->second;
            set<int>::iterator iter_set;
            int m=0;
            for(iter_set=s.begin();iter_set!=s.end();iter_set++){
                char *s = (char*)whole_gene+basic_start+(*iter_set);
                if(s==NULL) cout << "no" <<endl;
                string str(s);
                //cout << str << endl;
                r_set[m]=make_pair(*iter_set,str);//xä¸ºè¯¥åç¼€åœ¨è¯¥åˆ†ç»„ä¸­èµ·å§‹ä½ç½®
                m++;
            }
            sort(r_set,r_set+m,compare);
            for(int k=0;k<m;k++){
                int r = order_num[iter->first]+k;
                Ti[r] = r_set[k].first+1;//Tiæ˜¯è¯¥ç»„ä¸­æ’åºåçš„SAå€¼
                g[r_set[k].first+1] = order[r_set[k].first] + r + 1;//1æ˜¯åŠ ä¸Šè‡ªå·±
            }
        }


//        for(int x=1;x<=l;x++){
//            cout << g[x] << " ";
//        }
//        cout << endl;
        End = times(&tmsend);
        cout << "g time: " << (End - Start)/(double) clktck << endl;


        //è®¡ç®—æ–°çš„CSAçš„å€¼
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

        //æ›´æ–°depart_pos
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

=======
#include <iostream>
#include <windows.h>
#include <fstream>
#include <vector>
#include <cstdio>
#include <string>
#include <math.h>
#include <algorithm>
#include<ctime>


using namespace std;

bool compare(pair<int,string> p1,pair<int,string> p2);
int main()
{
    ifstream f;
    vector<char> whole_gene;
    string tmp;
    f.open("NC_008253.fna",ios::in);
    //f.open("testtext.fna",ios::in);

    if(f.fail()){
        cout << "can not open the file." << endl;
        exit(1);
    }
    //È¥µôµÚÒ»ĞĞ
    if(f.good())
        getline(f,tmp);
    do{
        if(f.good()){
            getline(f,tmp);
            for(string::iterator it=tmp.begin();it<tmp.end();it++){
                whole_gene.push_back(*it);
            }
        }
    }while(!f.eof());
    whole_gene.push_back('$');
    //cout << whole_gene.size() << endl;

    int n=whole_gene.size();//¼ÓÉÏ×îºóµÄ$
    int k=(int)5*log2(n);//×éÊı
    int l=n/k;//Ã¿×éµÄ×Ö·ûÊı
//    int l=6;
//    int k=4;

    int last_num = n-k*l;

    //¸ø×îºóÒ»×én-l(k-1) : nÅÅĞò²¢¼ÆËãCSA
    vector<pair<int,string> > T;
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
    cout << "T.size = " << T.size() << endl;


//    cout << "SA" << endl;
//    for(int i=0;i<tl;i++){
//        cout << T[i].first << " ";
//    }
//    cout << endl;
    //ÇóSA-1
    int RSA[tl];//SA-1
    for(int i=0;i<tl;i++){
        RSA[T[i].first] = i;
    }
//    cout << "SA-1" << endl;
//    for(int i=0;i<tl;i++){
//        cout << RSA[i] << " ";
//    }
//    cout << endl;
    //ÇóCSAÖµ
    vector<int> CSA(tl);
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

    //Çó·Ö½çµã
    int depart_pos[4] = {-1};//ÓÃ-1À´±ê¼Ç
    depart_pos[0] = 0;
    for(vector<pair<int,string>>::iterator iter=T.begin();iter!=T.end();iter++){
        if((*iter).second.at(0)=='A')
            depart_pos[1] = iter - T.begin();
        if((*iter).second.at(0)=='C')
            depart_pos[2] = iter - T.begin();
        if((*iter).second.at(0)=='G')
            depart_pos[3] = iter - T.begin();
    }


    //merge step
    int order[l+1];//×îºóÒ»¸ö´æµÄÊÇT¡®µÄorder
    order[l] = RSA[0];
    //cout << order[l] << endl;
    vector<pair<int,string> > Ti(l+1);
    for(int i=0;i<k;i++){//¶ÔÃ¿Ò»×é£¬´ÓºóÍùÇ°±àºÅ

//        cout << "depart_pos" << endl;
//        for(int i=0;i<4;i++)
//            cout << depart_pos[i] << " ";
//        cout << endl;

        int basic = n - i*l-last_num-1;//Ã¿Ò»×é×îºóÒ»¸ö×Ö·ûÔÚwhole_geneÖĞµÄÎ»ÖÃ
        cout << "handling "<<basic-l+1 <<" to "<<basic<<endl;
        int T_size = i*l + last_num;//ÓÃT±íÊ¾ÒÑ¾­ÅÅºÃĞòµÄ²¿·Ö
        int depart_pos_increase[4]={0};

        DWORD Start = timeGetTime();
        //¶Ô¸Ã×éÖĞµÄÃ¿Ò»¸öºó×º£¬´ÓºóÍùÇ°¼ÆËãorderÖµ
        for(int j=l;j>0;j--){
            char c = whole_gene[basic-l+j];
            switch(c){
                case 'A':
                    depart_pos_increase[1]++;
                    //depart_pos_increase[2]++;
                    //depart_pos_increase[3]++;
                    if(depart_pos[1]==-1 || CSA[depart_pos[0]+1]>order[j]){//B=¿Õ¼¯
                        order[j-1]= depart_pos[0];
                        break;
                    }
                    for(int k=depart_pos[0]+1;k<=depart_pos[1];k++){
                        if(CSA[k]<=order[j]){
                            order[j-1] = k;
                        }
                    }
                    break;
                case 'C':
                    depart_pos_increase[2]++;
                    //depart_pos_increase[3]++;
                    if(depart_pos[2]==-1 || CSA[depart_pos[1]+1]>order[j]){//B=¿Õ¼¯
                        order[j-1]= depart_pos[1];
                        break;
                    }
                    for(int k=depart_pos[1]+1;k<=depart_pos[2];k++){
                        if(CSA[k]<=order[j]){
                            order[j-1] = k;
                        }
                    }
                    break;
                case 'G':
                    depart_pos_increase[3]++;
                    if(depart_pos[3]==-1 || CSA[depart_pos[2]+1]>order[j]){//B=¿Õ¼¯
                        order[j-1]= depart_pos[2];
                        break;
                    }
                    for(int k=depart_pos[2]+1;k<=depart_pos[3];k++){
                        if(CSA[k]<=order[j]){
                            order[j-1] = k;
                        }
                    }
                    break;
                case 'T':
                    depart_pos_increase[0]++;
                    if(depart_pos[3]+1==T_size){//Ö®Ç°Ã»ÓĞÒÔT¿ªÍ·µÄºó×º
                        order[j-1]= depart_pos[3];
                        break;
                    }
                    if(CSA[depart_pos[3]+1]>order[j]){//B=¿Õ¼¯
                        order[j-1]= depart_pos[3];
                        break;
                    }
                    for(int k=depart_pos[3]+1;k<=i*l+last_num-1;k++){
                        if(CSA[k]<=order[j]){
                            order[j-1] = k;
                        }
                    }
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
        DWORD End = timeGetTime();
        cout << "order time: " << End-Start << endl;


        //½¨Á¢f(j)Ô­À´TÖĞÏÖÔÚµÄÖµ
        //±éÀúorder£¬¼ÆËãorderĞ¡ÓÚ¸÷¸öÖµµÄÊıÁ¿
        Start = timeGetTime();
        int order_num[T_size] = {0};
        for(int x=0;x<l;x++){
            for(int y=T_size-1;y>order[x];y--){
                order_num[y]++;
            }
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
        End = timeGetTime();
        cout << "f time: " << End-Start << endl;

        //½¨Á¢g(j)ĞÂ¼ÓµÄºó×ºµÄÖµ
        Start = timeGetTime();
        int g[l+1];
//        cout << "g" << endl;
        Ti.clear();
        int sa = l;//¸Ã×éÖĞ×îºóÒ»¸ö×Ö·ûµÄÎ»ÖÃ0-l
        for(int i=basic;i>basic-l;i--){
            s.insert(s.begin(),whole_gene[i]);
//            cout << s << endl;
            Ti.push_back(make_pair(sa,s));
            sa--;
        }
//        cout << Ti.capacity()<<endl;
        //cout << sa;
        sort(Ti.begin(),Ti.end(),compare);//ÅÅµÄÊÇÔ­À´±¾×éĞÂ¼ÓµÄl¸öºó×º
        for(int k=1;k<=l;k++){
            g[Ti[k-1].first] = order[Ti[k-1].first-1] + k;
        }
//        for(int x=1;x<=l;x++){
//            cout << g[x] << " ";
//        }
//        cout << endl;

//        for(int k=1;k<=l;k++){
//            //¼ÆËãµÚk¸öºó×ºÔÚ¸Ã×éÖĞµÄÅÅÃû
//            int r = order_num[order[k-1]];
//            for(int x=0;x<l;x++){
////                if(order[x]<order[k-1]){
////                    r++;
////                    continue;
////                }
//                if(order[x]==order[k-1]){
//                    //xºÍKÖ±½Ó±È¶Ô
//                    int minsize = l-max(x,k-1)+T_size;
//                    int x_start = basic-l+x+1;
//                    int k_start = basic-l+k;
//                    for(int y=0;y<minsize;y++){
//                        if(whole_gene[x_start+y]<whole_gene[k_start+y]){
//                            r++;
//                            break;
//                        }
//                        else if(whole_gene[x_start+y]>whole_gene[k_start+y])
//                            break;
//                    }
//                }
//            }
//            Ti[r]=make_pair(k,NULL);
//            g[k] = order[k-1] + r + 1;//1ÊÇ¼ÓÉÏ×Ô¼º
//        }

        End = timeGetTime();
        cout << "g time: " << End-Start << endl;


        //¼ÆËãĞÂµÄCSAµÄÖµ
        Start = timeGetTime();
        int T_size_new = T_size+l;
        vector<int> CSA_new(T_size_new);
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
                CSA_new[t] = g[Ti[jg-1].first+1];
                //cout << Ti[jg-1].first << " " << jg <<endl;
                jg++;
            }
        }
        End = timeGetTime();
        cout << "update CSA time: " << End-Start << endl;
        CSA = CSA_new;

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

        //¸üĞÂdepart_pos
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
    }
    for(int i=0;i<n;i++)
        cout << CSA[i] << " ";
    cout << endl;

    return 0;


}

bool compare(pair<int,string> p1, pair<int,string> p2){
    string s1 = p1.second;
    string s2 = p2.second;
	return s1+s2 < s2+s1;
}



>>>>>>> 8f33938aa99845ed0777145ae3c14ba80a7d0033
