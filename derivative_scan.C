#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <iomanip>

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"

using namespace std;

//============================================================
// Data structure
//============================================================

struct DataPoint{
    double ITHR;
    double IDB;
    double ICASN;
    double TH;
};

//============================================================
// String helpers
//============================================================

string trim(const string &s)
{
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end   = s.find_last_not_of(" \t\r\n");

    if(start==string::npos) return "";
    return s.substr(start,end-start+1);
}

bool safe_stod(const string& s,double& val)
{
    try{
        val = stod(trim(s));
        return true;
    }
    catch(...){
        return false;
    }
}

string format1(double val)
{
    stringstream ss;
    ss<<fixed<<setprecision(1)<<val;
    return ss.str();
}

//============================================================
// Derivative (central difference like numpy.gradient)
//============================================================

vector<double> computeDerivative(const vector<double>& x,
                                 const vector<double>& y)
{
    vector<double> dydx(x.size());

    for(size_t i=0;i<x.size();i++)
    {
        if(i==0)
            dydx[i]=(y[i+1]-y[i])/(x[i+1]-x[i]);

        else if(i==x.size()-1)
            dydx[i]=(y[i]-y[i-1])/(x[i]-x[i-1]);

        else
            dydx[i]=(y[i+1]-y[i-1])/(x[i+1]-x[i-1]);
    }

    return dydx;
}

//============================================================
// CSV Reader
//============================================================

vector<DataPoint> readCSV(string filename)
{
    vector<DataPoint> data;

    ifstream file(filename);

    if(!file.is_open())
    {
        cout<<"Cannot open "<<filename<<endl;
        return data;
    }

    string line;

    getline(file,line); // header

    while(getline(file,line))
    {
        stringstream ss(line);

        string s_ithr,s_idb,s_icasn,s_th,s_tmp;

        getline(ss,s_ithr,',');
        getline(ss,s_idb,',');
        getline(ss,s_icasn,',');
        getline(ss,s_th,','); // TH_MEAN

        double ithr,idb,icasn,th;

        if(!safe_stod(s_ithr,ithr)) continue;
        if(!safe_stod(s_idb,idb)) continue;
        if(!safe_stod(s_icasn,icasn)) continue;
        if(!safe_stod(s_th,th)) continue;

        data.push_back({ithr,idb,icasn,th});
    }

    cout<<"Loaded "<<data.size()<<" rows from "<<filename<<endl;

    return data;
}

//============================================================
// Plot helper
//============================================================

void makePlot(vector<double>& x_ir,
              vector<double>& y_ir,
              vector<double>& x_un,
              vector<double>& y_un,
              string xlabel,
              string ylabel,
              string title,
              string outname)
{

    TCanvas *c = new TCanvas(outname.c_str(),"",800,600);

    TGraph *g1 = new TGraph(x_ir.size(),&x_ir[0],&y_ir[0]);
    TGraph *g2 = new TGraph(x_un.size(),&x_un[0],&y_un[0]);

    //g1->SetMarkerStyle(20);
    //g1->SetLineWidth(2);

    //g2->SetMarkerStyle(24);
    //g2->SetLineWidth(2);

    // Irradiated (circle)
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(1.3);
    //g1->SetMarkerColor(kRed+1);
    //g1->SetLineColor(kRed+1);
    g1->SetMarkerColor(kOrange+7);
    g1->SetLineColor(kOrange+7);
    g1->SetLineWidth(2);

    // Unirradiated (triangle)
    g2->SetMarkerStyle(22);
    g2->SetMarkerSize(1.4);
    //g2->SetMarkerColor(kBlue+1);
    //g2->SetLineColor(kBlue+1);
    g2->SetMarkerColor(kAzure+2);
    g2->SetLineColor(kAzure+2);
    g2->SetLineWidth(2);

    g1->Draw("APL");
    g2->Draw("PL SAME");

    g1->GetXaxis()->SetTitle(xlabel.c_str());
    g1->GetYaxis()->SetTitle(ylabel.c_str());

    g1->SetTitle("");

    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.04);
    text.DrawLatex(0.15,0.94,title.c_str());

    TLegend *leg = new TLegend(0.65,0.8,0.88,0.88);
    leg->AddEntry(g1,"Irradiated","pl");
    leg->AddEntry(g2,"Unirradiated","pl");
    leg->Draw();

    c->SaveAs((outname+".pdf").c_str());
}

//============================================================
// Generic derivative scan
//============================================================

void derivativeScan(vector<DataPoint>& ir,
                    vector<DataPoint>& un,
                    string var,
                    string fixed1,
                    string fixed2)
{

    set<double> values1;
    set<double> values2;

    for(auto&p:ir)
    {
        if(fixed1=="ITHR") values1.insert(p.ITHR);
        if(fixed1=="IDB")  values1.insert(p.IDB);
        if(fixed1=="ICASN")values1.insert(p.ICASN);

        if(fixed2=="ITHR") values2.insert(p.ITHR);
        if(fixed2=="IDB")  values2.insert(p.IDB);
        if(fixed2=="ICASN")values2.insert(p.ICASN);
    }

    for(auto v1:values1)
    for(auto v2:values2)
    {

        vector<double> xir,yir,xun,yun;

        auto process=[&](vector<DataPoint>& data,
                         vector<double>& x,
                         vector<double>& y)
        {

            vector<pair<double,double>> points;

            for(auto&p:data)
            {

                double f1=(fixed1=="ITHR")?p.ITHR:
                          (fixed1=="IDB")?p.IDB:p.ICASN;

                double f2=(fixed2=="ITHR")?p.ITHR:
                          (fixed2=="IDB")?p.IDB:p.ICASN;

                if(f1!=v1 || f2!=v2) continue;

                double xv=(var=="ITHR")?p.ITHR:
                          (var=="IDB")?p.IDB:p.ICASN;

                points.push_back({xv,p.TH});
            }

            if(points.size()<2) return;

            sort(points.begin(),points.end());

            vector<double> xx,yy;

            for(auto&p:points)
            {
                xx.push_back(p.first);
                yy.push_back(p.second);
            }

            vector<double> dydx=computeDerivative(xx,yy);

            for(size_t i=0;i<xx.size();i++)
            {
                x.push_back(xx[i]);
                y.push_back(dydx[i]);
            }

        };

        process(ir,xir,yir);
        process(un,xun,yun);

        if(xir.size()<2) continue;

        string title=fixed1+"="+format1(v1)+"   "+fixed2+"="+format1(v2);

        string out="plots_derivatives/dTH_d"+var+"_"+fixed1+"_"+format1(v1)+"_"+fixed2+"_"+format1(v2);

        makePlot(xir,yir,xun,yun,
                 var,
                 "d(Threshold)/d"+var,
                 title,
                 out);
    }
}

//============================================================
// MAIN
//============================================================

void derivative_scan()
{

    vector<DataPoint> ir =
        readCSV("iitm_MALTA2testing-threshold_noise_irrediated.csv");

    vector<DataPoint> un =
        readCSV("iitm_MALTA2testing-threshold_noise_unirrediated.csv");

    derivativeScan(ir,un,"ITHR","IDB","ICASN");

    derivativeScan(ir,un,"IDB","ITHR","ICASN");

    derivativeScan(ir,un,"ICASN","IDB","ITHR");

}