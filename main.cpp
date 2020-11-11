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
    //ȥ����һ��
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

    int n=whole_gene.size();//��������$
    int k=(int)5*log2(n);//����
    int l=n/k;//ÿ����ַ���
//    int l=6;
//    int k=4;

    int last_num = n-k*l;

    //�����һ��n-l(k-1) : n���򲢼���CSA
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
    //��SA-1
    int RSA[tl];//SA-1
    for(int i=0;i<tl;i++){
        RSA[T[i].first] = i;
    }
//    cout << "SA-1" << endl;
//    for(int i=0;i<tl;i++){
//        cout << RSA[i] << " ";
//    }
//    cout << endl;
    //��CSAֵ
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

    //��ֽ��
    int depart_pos[4] = {-1};//��-1�����
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
    int order[l+1];//���һ�������T����order
    order[l] = RSA[0];
    //cout << order[l] << endl;
    vector<pair<int,string> > Ti(l+1);
    for(int i=0;i<k;i++){//��ÿһ�飬�Ӻ���ǰ���

//        cout << "depart_pos" << endl;
//        for(int i=0;i<4;i++)
//            cout << depart_pos[i] << " ";
//        cout << endl;

        int basic = n - i*l-last_num-1;//ÿһ�����һ���ַ���whole_gene�е�λ��
        cout << "handling "<<basic-l+1 <<" to "<<basic<<endl;
        int T_size = i*l + last_num;//��T��ʾ�Ѿ��ź���Ĳ���
        int depart_pos_increase[4]={0};

        DWORD Start = timeGetTime();
        //�Ը����е�ÿһ����׺���Ӻ���ǰ����orderֵ
        for(int j=l;j>0;j--){
            char c = whole_gene[basic-l+j];
            switch(c){
                case 'A':
                    depart_pos_increase[1]++;
                    //depart_pos_increase[2]++;
                    //depart_pos_increase[3]++;
                    if(depart_pos[1]==-1 || CSA[depart_pos[0]+1]>order[j]){//B=�ռ�
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
                    if(depart_pos[2]==-1 || CSA[depart_pos[1]+1]>order[j]){//B=�ռ�
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
                    if(depart_pos[3]==-1 || CSA[depart_pos[2]+1]>order[j]){//B=�ռ�
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
                    if(depart_pos[3]+1==T_size){//֮ǰû����T��ͷ�ĺ�׺
                        order[j-1]= depart_pos[3];
                        break;
                    }
                    if(CSA[depart_pos[3]+1]>order[j]){//B=�ռ�
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


        //����f(j)ԭ��T�����ڵ�ֵ
        //����order������orderС�ڸ���ֵ������
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

        //����g(j)�¼ӵĺ�׺��ֵ
        Start = timeGetTime();
        int g[l+1];
////        cout << "g" << endl;
//        Ti.clear();
//        int sa = l;//���������һ���ַ���λ��0-l
//        for(int i=basic;i>basic-l;i--){
//            s.insert(s.begin(),whole_gene[i]);
////            cout << s << endl;
//            Ti.push_back(make_pair(sa,s));
//            sa--;
//        }
////        cout << Ti.capacity()<<endl;
//        //cout << sa;
//        sort(Ti.begin(),Ti.end(),compare);//�ŵ���ԭ�������¼ӵ�l����׺
//        for(int k=1;k<=l;k++){
//            g[Ti[k-1].first] = order[Ti[k-1].first-1] + k;
//        }
////        for(int x=1;x<=l;x++){
////            cout << g[x] << " ";
////        }
////        cout << endl;

        for(int k=1;k<=l;k++){
            //�����k����׺�ڸ����е�����
            int r = order_num[order[k-1]];
            for(int x=0;x<l;x++){
//                if(order[x]<order[k-1]){
//                    r++;
//                    continue;
//                }
                if(order[x]==order[k-1]){
                    //x��Kֱ�ӱȶ�
                    int minsize = l-max(x,k-1)+T_size;
                    int x_start = basic-l+x+1;
                    int k_start = basic-l+k;
                    for(int y=0;y<minsize;y++){
                        if(whole_gene[x_start+y]<whole_gene[k_start+y]){
                            r++;
                            break;
                        }
                        else if(whole_gene[x_start+y]>whole_gene[k_start+y])
                            break;
                    }
                }
            }
            Ti[r]=make_pair(k,NULL);
            g[k] = order[k-1] + r + 1;//1�Ǽ����Լ�
        }

        End = timeGetTime();
        cout << "g time: " << End-Start << endl;


        //�����µ�CSA��ֵ
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

        //����depart_pos
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



