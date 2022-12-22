#include "../include/tskf_estim.h"

void TSKF::load_parameters() {

    if( !_nh.getParam("motor_num", _motor_num) ) {
        _motor_num =  4;
    }

    if( !_nh.getParam("model_name", _model_name) ) {
        _model_name =  "iris_smc";
    }
    vector<double> inertia;
    if( !_nh.getParam("inertia", inertia) ) {
        inertia.resize(3);
        inertia[0] = 1.0;
        inertia[1] = 1.0;
        inertia[2] = 1.0;
    }
    _inertia = Eigen::Matrix3d( Eigen::Vector3d( inertia[0], inertia[1], inertia[2] ).asDiagonal() );

    if( !_nh.getParam( "arm_length", _arm_length ) ) {
        _arm_length.resize( _motor_num );
        _arm_length[0] = 0.255;
        _arm_length[1] = 0.238;
        _arm_length[2] = 0.255;
        _arm_length[3] = 0.238;
    }

    if( !_nh.getParam( "motor_rotation_direction", _motor_rotation_direction ) ) {
        _motor_rotation_direction.resize( _motor_num );
        _motor_rotation_direction[0] = 1;
        _motor_rotation_direction[1] = 1;
        _motor_rotation_direction[2] = -1;
        _motor_rotation_direction[3] = -1;
    }

    if( !_nh.getParam("motor_force_k", _motor_force_k) ) {
        _motor_force_k =  8.54858e-06;
    }
    if( !_nh.getParam("motor_moment_k", _motor_moment_k) ) {
        _motor_moment_k = 1.6e-2;
    }

    if( !_nh.getParam("mass", _mass) ) {
        _mass =  1.5;
    }

    if( !_nh.getParam("gravity", _gravity) ) {
        _gravity =  9.81;
    }

    if( !_nh.getParam( "rotor_angles", _rotor_angles ) ) {
        _rotor_angles.resize( _motor_num );
        _rotor_angles[0] = -0.5337;
        _rotor_angles[1] = 2.565;
        _rotor_angles[2] = 0.5337;
        _rotor_angles[3] = -2.565;
    }

}

TSKF::TSKF() {

    load_parameters();

    //Initialization vectors and matrices
    _x_hat << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
      
    _x_tilde << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    _gamma << 0.0, 0.0, 0.0, 0.0;

    _res = Eigen::Matrix<double,6,1>::Zero();

    _V_kk_kk = Eigen::Matrix<double,12,4>::Zero();
    _P_gamma_kk_kk = Eigen::Matrix<double,4,4>::Zero();
    _P_x_kk_kk = Eigen::Matrix<double,12,12>::Zero();

    _vec << 0.0, 0.0, 0.0, 0.0, 0.0, -9.81*_Ts, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0;

    _C_k << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    _A_k << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, -_gravity, 0, 0, 0, 0, 0,
     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, -_gravity, 0, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    _odom_sub = _nh.subscribe(_model_name + "/odometry", 0, &TSKF::odometry_cb, this); 
    _mot_vel_sub = _nh.subscribe < std_msgs::Float32MultiArray > (_model_name + "/cmd/motor_vel", 0, &TSKF::mot_vel_cb, this);
    _motori = _nh.advertise<std_msgs::Float32MultiArray>("/cmd/motor_filter", 0);   
    _yy = _nh.advertise<std_msgs::Float32MultiArray>( "/y_kk", 0 );   
}

void TSKF::odometry_cb (const nav_msgs::Odometry odometry_msg) {
    _odom = odometry_msg;
    Eigen::Vector3d p;
    Eigen::Vector3d w;

    p << _odom.pose.pose.position.x, _odom.pose.pose.position.y, _odom.pose.pose.position.z;
    w << _odom.twist.twist.angular.x, _odom.twist.twist.angular.y, _odom.twist.twist.angular.z;

    change_y( p, w );
}

void TSKF::change_y( Eigen::Vector3d p, Eigen::Vector3d w ) { //corretto
    Eigen::Vector3d p_new;
    Eigen::Vector3d w_new;

    p_new = _Rx*p;
    w_new = _Rx*w;

    _y_kk << p_new[0], p_new[1], p_new[2], w_new[1], w_new[0], w_new[2];
}

void TSKF::mot_vel_cb (const std_msgs::Float32MultiArray mot_vel) {
    _mot = mot_vel;
    Eigen::Vector4d motor_vel;

    motor_vel << _mot.data[0], _mot.data[1], _mot.data[2], _mot.data[3];
    change_u( motor_vel );
}

void TSKF::change_u ( Eigen::Vector4d motor_vel ) {
    _u_k << motor_vel[0], motor_vel[1], motor_vel[2], motor_vel[3];
    _u_k = _u_k*0.01;
}


void TSKF::publisher_test() {

    ros::Rate r(100);
    std_msgs::Float32MultiArray veloc;
    std_msgs::Float32MultiArray ykk;
    veloc.data.resize( _motor_num );
    ykk.data.resize(6);

    while( ros::ok() ) {

        for( int i=0; i<_motor_num; i++ ) {
            veloc.data[i] = _u_k(i);
        }

        for( int i=0; i<6; i++ ) {
            ykk.data[i] = _y_kk(i);
        }

        _motori.publish( veloc );
        _yy.publish( ykk );

        r.sleep();
    }
}

bool TSKF::generate_allocation_matrix(Eigen::MatrixXd & allocation_M, 
                                    int motor_size,
                                    vector<double> rotor_angle,
                                    vector<double> arm_length, 
                                    double force_k,
                                    double moment_k,
                                    vector<int> direction ) {

    allocation_M.resize(4, motor_size );

    for(int i=0; i<motor_size; i++ ) {
        allocation_M(0, i) = sin( rotor_angle[i] ) * arm_length[i] * force_k;
        allocation_M(1, i) = cos( rotor_angle[i] ) * arm_length[i] * force_k;
        allocation_M(2, i) = direction[i] * force_k * moment_k;
        allocation_M(3, i) = -force_k;
    }

    Eigen::FullPivLU<Eigen::Matrix4Xd> lu( allocation_M);
    if ( lu.rank() < 4 ) {
        ROS_ERROR("The allocation matrix rank is lower than 4. This matrix specifies a not fully controllable system, check your configuration");
        return false;
    }

    return true;
}

void TSKF::tskf_matrix_generation( Eigen::MatrixXd allocation_matrix ) {
    Eigen::Matrix<double,12,12> I_12;           //matrice identità
    _A_k = I_12 + _A_k*_Ts;
}

void TSKF::estimation(Eigen::Vector4d u, Eigen::MatrixXd y, Eigen::Vector4d pr_gamma) {
    
    ros::Rate r( 100 );

    Eigen::MatrixXd allocation_M;
    if(!generate_allocation_matrix( allocation_M, _motor_num, _rotor_angles, _arm_length, _motor_force_k, _motor_moment_k, _motor_rotation_direction ) ) {     
        cout << "Wrong allocation matrix" << endl;
        exit(0);
    }
    //Variabili del tskf usate solo in estimation
    Eigen::Matrix<double,12,4> W_k;             //coupling
    Eigen::Matrix<double,4,4> P_gamma_kk_k;     //predizione P_gamma
    Eigen::Vector4d gamma_pred;                 //predizione gamma
    Eigen::Matrix<double,12,4> V_kk_k;          //coupling
    Eigen::Matrix<double,6,4> H_kk_k;           //coupling
    Eigen::Matrix<double,12,12> P_x_kk_k;       //predizione P_x
    Eigen::Matrix<double,12,1> x_kk_k;          //predizione x
    Eigen::Matrix<double,12,6> Kx_kk;           //guadagno kf per lo stato
    Eigen::Matrix<double,6,6> S_kk;             //residuo
    Eigen::Matrix<double,4,6> K_gamma_kk;       //guadagno kf per gamma
    Eigen::Matrix<double,4,4> U;                //input di controllo (matrice)
    Eigen::Matrix<double,4,4> I_4;              //matrice identità
    Eigen::Matrix<double,12,12> I_12;           //matrice identità

    I_4.diagonal() << 1.0, 1.0, 1.0, 1.0;
    I_12.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    
    U << u[0], 0.0, 0.0, 0.0,
    0.0, u[1], 0.0, 0.0,
    0.0, 0.0, u[2], 0.0,
    0.0, 0.0, 0.0, u[3];

    while (ros::ok() ) {

        W_k = _A_k*_V_kk_kk - _B_k*U; //coupling equation (26)
        P_gamma_kk_k = _P_gamma_kk_kk + _Qgamma; //optimal bias estimator (15)
        gamma_pred = _gamma; //OBE (14)
        V_kk_k = W_k*_P_gamma_kk_kk*P_gamma_kk_k.inverse(); //coupling (27)
        H_kk_k = _C_k*V_kk_k; //coupling (28)
        P_x_kk_k = _A_k*_P_x_kk_kk*_A_k.transpose() + _Qx + W_k*_P_gamma_kk_kk*W_k.transpose()
        - V_kk_k*P_gamma_kk_k*V_kk_k.transpose(); //BFE (20)
        x_kk_k = _A_k*_x_tilde + W_k*_gamma - V_kk_k*_gamma + _vec; //BFE (19)
        _res = y - _C_k*x_kk_k; //calcolo residui (24)
        Kx_kk = P_x_kk_k*_C_k.transpose()*(_C_k*P_x_kk_k*_C_k.transpose() + _R).inverse(); //BFE (22)
        S_kk = _C_k*P_x_kk_k*_C_k.transpose() + _R; //calcolo residui (25)
        K_gamma_kk = P_gamma_kk_k*(H_kk_k.transpose())*(H_kk_k*P_gamma_kk_k*H_kk_k.transpose() + S_kk).inverse(); //OBE (17)

        //output stimatore
        _V_kk_kk = V_kk_k - Kx_kk*H_kk_k; //coupling (29)
        _P_gamma_kk_kk = (I_4 - K_gamma_kk*H_kk_k)*P_gamma_kk_k; //OBE (18)
        _P_x_kk_kk = (I_12 - Kx_kk*_C_k)*P_x_kk_k; //BFE (23)
        _x_tilde = x_kk_k + Kx_kk*( y - _C_k*x_kk_k ); //BFE (21)
        _gamma = gamma_pred + K_gamma_kk*( _res - H_kk_k*pr_gamma); //OBE (16)

        r.sleep();

    }


}

void TSKF::run() {

    boost::thread publishr_test_t( &TSKF::publisher_test, this );
    ros::spin();
}

int main( int argc, char** argv ) {

    ros::init(argc, argv, "tskf_estimator");
    TSKF e;
    e.run();

    return 0;
}