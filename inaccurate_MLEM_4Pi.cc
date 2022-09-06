#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <list>
#include <vector>
#include <time.h>
#include <assert.h>

using namespace std;

#define Electron_Mass 510.99
#define GammaEnergy 511.0 //keV
#define Electron_Radius 2.818 //fm
#define Pi 3.1415
#define Rec_Theta_Bins 180
#define Rec_Phi_Bins 360
#define Cal_Theta_Bins 180
#define Cal_Phi_Bins 360
#define Bot_Theta 0.0
#define Top_Theta 180.0
#define Bot_Phi 0.0
#define Top_Phi 360.0
#define MLEM 2


class Density
{
  public:
  Density(int theta_num, int phi_num); 
  Density();
  ~Density() {};
  double output_count(int theta_index, int phi_index);
  void input_count(int theta_index, int phi_index);
  void input_system(int theta_index, int phi_index, double den);
  void change(int theta_index, int phi_index, double den);

  public:
  double** count;
  int theta_number, phi_number;
};
Density::Density()
{}
Density::Density(int theta_num, int phi_num)
{
  count = new double*[theta_num];
  for (int i = 0; i < theta_num; i++)
  {
    count[i]=new double[phi_num];
  }

  for (int i = 0; i < theta_num; i++)
  for (int j = 0; j < phi_num; j++)
  {
    count[i][j]=0;
  }

  theta_number = theta_num;
  phi_number = phi_num;
}
void Density::input_count(int theta_index, int phi_index)
{
  if ((theta_index < theta_number) and (phi_index < phi_number))
  {
    count[theta_index][phi_index]++;
  }
}
void Density::input_system(int theta_index, int phi_index, double den)
{
  if ((theta_index < theta_number) and (phi_index < phi_number))
  {
    count[theta_index][phi_index] += den;
  }
}
void Density::change(int theta_index, int phi_index, double den)
{
  if ((theta_index < theta_number) and (phi_index < phi_number))
  {
    count[theta_index][phi_index] = den;
  }
}
double Density::output_count(int theta_index, int phi_index)
{
  return count[theta_index][phi_index];
}

// define Compton event class
class Compton_event
{
public:
  Compton_event(double Com_axes_theta, double Com_axes_phi, double Com_scatter_theta); 
  Compton_event();
  ~Compton_event() {};
  double ComptonCrossSection();

public:
  double Com_ax_theta, Com_ax_phi, Com_sca_theta;
  double cos_Com_sca_theta;
  double sin_Com_sca_theta;
  double cos_Com_ax_theta;
  double sin_Com_ax_theta;
};
Compton_event::Compton_event()
{}
Compton_event::Compton_event(double Com_axes_theta, double Com_axes_phi, double Com_scatter_theta)
{
  Com_ax_theta = Com_axes_theta;
  Com_ax_phi = Com_axes_phi;
  Com_sca_theta = Com_scatter_theta;

  cos_Com_sca_theta = cos((Com_sca_theta/180)*Pi);
  sin_Com_sca_theta = sin((Com_sca_theta/180)*Pi);
  cos_Com_ax_theta = cos((Com_ax_theta/180)*Pi);
  sin_Com_ax_theta = sin((Com_ax_theta/180)*Pi);
}
double Compton_event::ComptonCrossSection()
{
  /***********Use the Klein-Nishina differential section formula***********/
  float Cos_Angle1 = cos_Com_sca_theta;
  float arfc = GammaEnergy/Electron_Mass;
  float R_gamma = 1.0/(1 + arfc*(1 - Cos_Angle1));
  float CrossSection1 = 0.5*Electron_Radius*Electron_Radius*R_gamma;
  float CrossSection2 = R_gamma*R_gamma + (Cos_Angle1*Cos_Angle1 - 1)*R_gamma + 1;
  float d_CrossSection = CrossSection1*CrossSection2;
  float TotalCrossSection = 2*Pi*sin_Com_sca_theta*d_CrossSection;
  return TotalCrossSection;
}

// simple check
bool Check_Sign(double a, double b, double c, double d)
{
  if ((a == 0) || (b == 0) || (c == 0) || (d == 0))
  {
    return true;
  }
  else
  {
    if ((a*b < 0) || (a*c < 0) || (a*d < 0) || (b*c < 0) || (b*d < 0) || (c*d < 0))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
}

bool Check_Degree(double theta, double phi, double cos_thetac, double cos_theta0, double sin_theta0, double phi_0, double phi_voxel, double theta_voxel) 
{
  double phi_top = phi + phi_voxel/2.0;
  double phi_bot = phi - phi_voxel/2.0;
  double theta_top = theta + theta_voxel/2.0;
  double theta_bot = theta - theta_voxel/2.0;
  
  double cos_theta_top = (cos_thetac - cos_theta0*cos((theta_top/180.0)*Pi))/(sin_theta0*sin((theta_top/180.0)*Pi));
  double cos_theta_bot = (cos_thetac - cos_theta0*cos((theta_bot/180.0)*Pi))/(sin_theta0*sin((theta_bot/180.0)*Pi));

  if (abs(cos_theta_top) > 1.00001 and abs(cos_theta_bot) > 1.00001) return false;
  double acos_theta_top = (acos(fminl(fmaxl(cos_theta_top,-1.0),1.0)))/Pi*180;
  double acos_theta_bot = (acos(fminl(fmaxl(cos_theta_bot,-1.0),1.0)))/Pi*180;

  // positive direction
  // phi = phi_0 + acos(right)
  double pos_phi_top = phi_0 + acos_theta_top; 
  double pos_phi_bot = phi_0 + acos_theta_bot; 

  double positive_left_top = phi_bot - pos_phi_top;
  double positive_right_top = phi_top - pos_phi_top;
  double positive_right_bot = phi_top - pos_phi_bot;
  double positive_left_bot = phi_bot - pos_phi_bot;
  bool positive_check = Check_Sign(positive_left_top, positive_right_top, positive_right_bot, positive_left_bot);

  // negative direction
  // phi = phi_0 - acos(right)
  double neg_phi_top = phi_0 - acos_theta_top;
  double neg_phi_bot = phi_0 - acos_theta_bot;

  double negative_left_top = phi_bot - neg_phi_top;
  double negative_right_top = phi_top - neg_phi_top;
  double negative_right_bot = phi_top - neg_phi_bot;
  double negative_left_bot = phi_bot - neg_phi_bot;
  bool negative_check = Check_Sign(negative_left_top, negative_right_top, negative_right_bot, negative_left_bot);

  if (positive_check || negative_check)
  {
    return true;
  }
  else
  {
    return false;
  }
}

int main()
{
    // read file
    ifstream infile;
    infile.open("Four_pi_data.txt",ios::in);
    cout<<" ifstream good!"<<endl;
    
    clock_t t_start;
    int rec_theta_bins =  Rec_Theta_Bins;
    int rec_phi_bins =  Rec_Phi_Bins;
    int cal_theta_bins =  Cal_Theta_Bins;
    int cal_phi_bins =  Cal_Phi_Bins;

    double cal_theta_voxel = (Top_Theta - Bot_Theta)/Cal_Theta_Bins;
    double cal_phi_voxel = (Top_Phi - Bot_Phi)/Cal_Phi_Bins;
    int factor_1 = cal_theta_bins/rec_theta_bins;
    int factor_2 = cal_phi_bins/rec_phi_bins;
    if (factor_1 != factor_2)
    {
      cout<<"error"<<endl;
      exit(1);
    }
    int factor = factor_1;

    double rec_theta_voxel = (Top_Theta - Bot_Theta)/Rec_Theta_Bins;
    double rec_phi_voxel = (Top_Phi - Bot_Phi)/Rec_Phi_Bins;

    double Com_axes_theta;
    double Com_axes_phi;
    double Com_scatter_theta;
    int count_event = 0;

    std::list<Compton_event*> Compton_event_Set;
    while (!infile.eof())
    {
        infile>>Com_axes_theta>>Com_axes_phi>>Com_scatter_theta;

        Compton_event* event= new Compton_event(Com_axes_theta, Com_axes_phi, Com_scatter_theta);
        Compton_event_Set.push_back(event);
        count_event++;
        // if (count_event > 1000) break;
    }
    std::list<Compton_event*>::iterator it1 = Compton_event_Set.end();
    it1--;
    it1 = Compton_event_Set.erase (it1);

    //construct count space and system matrix
    Density recon_space(rec_theta_bins, rec_phi_bins);
    Density SystemMatrix(rec_theta_bins, rec_phi_bins);

    // print out data
    int len = Compton_event_Set.size();
    cout<<"Compton Event number: "<<len<<endl;

    // fill the density matrix
    std::list<Compton_event*>::iterator itr;
    double cos_thetac, cos_theta0, sin_theta0;
    double phi_0, theta_0, thetac;
    for(itr=Compton_event_Set.begin(); itr!=Compton_event_Set.end(); itr++)
    {
      if ((*itr)->Com_ax_theta == 0)
      {
        double re_theta = (*itr)->Com_sca_theta;
        int theta_index = int(re_theta/rec_theta_voxel);
        for (int i = 0; i < rec_phi_bins; i++)
        {
          recon_space.input_count(theta_index, i);
          SystemMatrix.input_system(theta_index, i, ((*itr)->ComptonCrossSection()));
        }
      }
      else if((*itr)->Com_ax_theta == 180)
      {
        double re_theta = 180 - (*itr)->Com_sca_theta;
        int theta_index = int(re_theta/rec_theta_voxel);
        for (int i = 0; i < rec_phi_bins; i++)
        {
          recon_space.input_count(theta_index,i);
          SystemMatrix.input_system(theta_index, i, ((*itr)->ComptonCrossSection()));
        }
      }
      else
      {
        cos_thetac = (*itr)->cos_Com_sca_theta;
        cos_theta0 = (*itr)->cos_Com_ax_theta;
        sin_theta0 = (*itr)->sin_Com_ax_theta;

        phi_0 = (*itr)->Com_ax_phi;
        theta_0 = (*itr)->Com_ax_theta;
        thetac = (*itr)->Com_sca_theta;

        int ini_theta = int(theta_0 - thetac) - 2;
        int end_theta = int(theta_0 + thetac) + 2;
        int ini_phi = int(phi_0 - 180);
        int end_phi = int(phi_0 + 180);
        int j_j;

        if (ini_theta < 0) ini_theta = 0;
        if (end_theta > 180) end_theta = 180;

        // projection for interested region
        for (int i = ini_theta; i < end_theta; i++)
        {
          for (int j = ini_phi; j < end_phi; j++)
          {
            
            double theta = i*cal_theta_voxel + cal_theta_voxel/2.0;
            double phi = j*cal_phi_voxel + cal_phi_voxel/2.0;
            bool judge = Check_Degree(theta, phi, cos_thetac, cos_theta0, sin_theta0, phi_0, cal_phi_voxel, cal_theta_voxel);
            if (judge)
            {
              j_j = j;
              int rec_theta_index = i/factor;
              if (j < 0) j += 360;
              if (j > 359) j -= 360;
              int rec_phi_index = j/factor;
              recon_space.input_count(rec_theta_index, rec_phi_index);
              SystemMatrix.input_system(rec_theta_index, rec_phi_index, ((*itr)->ComptonCrossSection()));
              j = j_j;
            }
          }
        }
      }
    }
    // cout<<"Simple back projection:  "<<float((clock() - t_start))/CLOCKS_PER_SEC<<"s"<<endl;

        // write out into file
    std::ofstream outList_1;
    outList_1.open("Simple_back_projection.txt",ios::out);
    if(!outList_1)
    {
      cout<<"Open the file failure..\n"<<endl;
    }
    for(int i=0;i<rec_theta_bins;i++)
    {
      for(int j=0;j<rec_phi_bins;j++)
      {
        outList_1<<recon_space.output_count(i,j)<<" ";
      }
      outList_1<<endl;
    }
    outList_1.close();  

    /*********************Normalization of SystemMatrix*********************/
    double SystemMatrixSum = 0;
    for(int Istheta = 0; Istheta < rec_theta_bins; Istheta++)
    {
      for(int Isphi = 0; Isphi < rec_phi_bins; Isphi++)
      {
        SystemMatrixSum += SystemMatrix.output_count(Istheta, Isphi);
      }
    }
    for(int Istheta = 0; Istheta < rec_theta_bins; Istheta++)
    {
      for(int Isphi = 0; Isphi < rec_phi_bins; Isphi++)
      {
          SystemMatrix.change(Istheta, Isphi,((SystemMatrix.output_count(Istheta, Isphi))/SystemMatrixSum));
      }
    }

    /************************MLEM itreation**************************************/
    for(int k = 0; k < MLEM; k++)
    {
      double TDetEvent=0;
      double FDetEvent=0;
      for(int Istheta = 0; Istheta < rec_theta_bins; Istheta++)
      {
        for(int Isphi = 0; Isphi < rec_phi_bins; Isphi++)
        {
          TDetEvent += recon_space.output_count(Istheta, Isphi)*SystemMatrix.output_count(Istheta, Isphi);
          FDetEvent += recon_space.output_count(Istheta, Isphi);
        }
      }
    
      Density CFactor(rec_theta_bins, rec_phi_bins);
      for(int Istheta = 0; Istheta < rec_theta_bins; Istheta++)
      {
        for(int Isphi = 0; Isphi < rec_phi_bins; Isphi++)
        {
          CFactor.change(Istheta, Isphi, (SystemMatrix.output_count(Istheta, Isphi)/TDetEvent));
          recon_space.change(Istheta, Isphi, (recon_space.output_count(Istheta, Isphi)*FDetEvent*CFactor.output_count(Istheta, Isphi)));
        } 
      }
   }
    // cout<<"MLEM iteration:  "<<float((clock() - t_start))/CLOCKS_PER_SEC<<"s"<<endl;

    // write out into file
    std::ofstream outList;
    outList.open("MLEM_back_projection.txt",ios::out);
    if(!outList)
    {
      cout<<"Open the file failure..\n"<<endl;
    }
    for(int i=0;i<rec_theta_bins;i++)
    {
      for(int j=0;j<rec_phi_bins;j++)
      {
        outList<<recon_space.output_count(i,j)<<" ";
      }
      outList<<endl;
    }
    outList.close();   
}