#include "ros/ros.h"
#include "nav_msgs/Odometry.h"
#include "boost/thread.hpp"
#include <Eigen/Eigen>
#include "gazebo_msgs/ModelStates.h"
#include "std_msgs/Float32MultiArray.h"
#include "std_msgs/Bool.h"
#include "tf/tf.h"
#include <Eigen/Geometry>
#include <fstream>
#include <iostream>
#include "utils.h"

using namespace std;
using namespace Eigen;

class TSKF {

    public:

        TSKF();

        void set_matrix(double mass_, Eigen::Matrix3d J_, double gravity_, double K_);
        //void estimation(Eigen::Vector4d pwm, Eigen::MatrixXd y, Eigen::Vector4d pr_gamma);
        void estimation();
        void change_y( Eigen::Vector3d p, Eigen::Vector3d w );
        void change_u( Eigen::Vector4d motor_vel );
        void odometry_cb(const nav_msgs::Odometry odometry_msg);
        void mot_vel_cb(const std_msgs::Float32MultiArray mot_vel);
        void publisher_test();
        bool generate_allocation_matrix(Eigen::MatrixXd & allocation_M, 
                                    int motor_size,
                                    vector<double> rotor_angle,
                                    vector<double> arm_length, 
                                    double force_k,
                                    double moment_k,
                                    vector<int> direction );
        void tskf_matrix_generation( Eigen::MatrixXd allocation_matrix); 
        void pwm_computation( );
        void controller_cb( const std_msgs::Float32MultiArray t_a_msg );
        void load_parameters();
        void estimator_reset();
        void reset_req_cb( std_msgs::Bool );
        void fault_generated_cb( const std_msgs::Float32MultiArray );
        void write_file();
        void run();

    private:
        
        ros::NodeHandle _nh;
        ros::Subscriber _odom_sub;
        ros::Subscriber _mot_vel_sub;
        ros::Subscriber _t_acc_sub;
        ros::Subscriber _estimator_reset_req;
        ros::Subscriber _f_generated;
        ros::Publisher _estimator_active;
        string _model_name;
        nav_msgs::Odometry _odom;
        std_msgs::Float32MultiArray _mot;
        vector<double> _arm_length;
        vector<int> _motor_rotation_direction;
        int _motor_num;
        double _motor_force_k;
        double _motor_moment_k;
        double _mass;
        double _gravity;
        Eigen::Matrix3d _inertia;
        vector<double> _rotor_angles;

        //publisher per debug
        ros::Publisher _motori;
        ros::Publisher _yy;
        ros::Publisher _x_estim;
        ros::Publisher _fault_estim;
        ros::Publisher _residual;
        ros::Publisher _fault_detection;

        //Estimator Inputs
        Eigen::Vector4d _u_k;
        Eigen::Matrix<double,6,1> _y_kk;
        
        //Estimator Outputs
        Eigen::Matrix<double,12,1>  _x_hat;
        //in the loop
        Eigen::Matrix<double,12,1> _x_tilde;
        Eigen::Vector4d _gamma;

        //Matrices used in TSKF

        //in the loop
        Eigen::Matrix<double,12,4> _V_kk_kk;
        Eigen::Matrix<double,4,4> _P_gamma_kk_kk;
        Eigen::Matrix<double,12,12> _P_x_kk_kk;
        Eigen::Matrix<double,6,1> _res;

        Eigen::Matrix<double,12,1> _vec;

        //Matrices modello linearizzato
        Eigen::Matrix<double,12,12> _A_k;
        Eigen::Matrix<double,12,4> _B_k;
        Eigen::Matrix<double,6,12> _C_k;

        //Matrici probabilistiche
        Eigen::Matrix<double,12,12> _Qx;
        Eigen::Matrix<double,4,4> _Qgamma;
        Eigen::Matrix<double,6,6> _R;
        
        Eigen::Matrix<double,4,4> _Ap_inv; //allocation matrix from the paper
        Eigen::Vector4d _u_p;
        Eigen::Vector4d _t_a;
        Eigen::Matrix4d _G;
        Eigen::Matrix4d _t2pwm;
        Eigen::Matrix<double,4,1> _detection;
        Eigen::Matrix<double,4,1> _generation;

        // Eigen::Vector4d _wd2rpm_new = Eigen::Matrix<double,4,1>::Zero();

        //Matrice cambio riferimento
        Eigen::Matrix3d _Rx = utilities::rotx(3.14);

        double _estim_rate;
        double _Ts;
        double _par;

        double _K = 175.0;
        double _Kpsi = 0.023;
        double _L = 0.2;

        double _r1;
        double _r2;
        double _qx1;
        double _qx2;
        double _qgamma;
        bool _est_res = false;
        



        




};