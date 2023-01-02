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
    _Qx = Eigen::Matrix<double,12,12>::Zero();
    _Qgamma = Eigen::Matrix<double,4,4>::Zero();
    _R = Eigen::Matrix<double,6,6>::Zero();

    _Qgamma.diagonal() << 0.1, 0.1, 0.1, 0.1;

    _Qx.diagonal() << 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6;

    _R.diagonal() << 1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6;

    _vec << 0.0, 0.0, 0.0, 0.0, 0.0, -9.81*_Ts, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0;

    _C_k << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    _A_k << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, _gravity, 0, 0, 0, 0, 0,
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

    _B_k << 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, //fino a qui non cambia nulla
    0.0, 0.0, 0.0, 0.0, //qui ci va il corrispettivo di K/m
    0.0, 0.0, 0.0, 0.0, //restano tutti zeri
    0.0, 0.0, 0.0, 0.0, //vanno modificate la colonna 1 e 3
    0.0, 0.0, 0.0, 0.0, //restano tutti zeri
    0.0, 0.0, 0.0, 0.0, //vanno modificati colonna 2 e 4
    0.0, 0.0, 0.0, 0.0, //restano tutti zeri
    0.0, 0.0, 0.0, 0.0, //vengono modificate tutte le colonne ma vedere bene motori!!

    _odom_sub = _nh.subscribe(_model_name + "/odometry", 0, &TSKF::odometry_cb, this); 
    _mot_vel_sub = _nh.subscribe < std_msgs::Float32MultiArray > (_model_name + "/cmd/motor_vel", 0, &TSKF::mot_vel_cb, this);
    _motori = _nh.advertise<std_msgs::Float32MultiArray>("/cmd/motor_filter", 0);   
    _yy = _nh.advertise<std_msgs::Float32MultiArray>( "/y_kk", 0 );   
    _x_estim = _nh.advertise<std_msgs::Float32MultiArray>( "/estim_states", 0 );
    _fault_estim = _nh.advertise<std_msgs::Float32MultiArray>( "/fault_estim", 0 );
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
    //_u_k = _u_k*_par;
}


void TSKF::publisher_test() {

    ros::Rate r(100);
    std_msgs::Float32MultiArray veloc;
    std_msgs::Float32MultiArray ykk;
    std_msgs::Float32MultiArray x_stim;
    std_msgs::Float32MultiArray f_stim;
    x_stim.data.resize(12);
    veloc.data.resize( _motor_num );
    ykk.data.resize(6);
    f_stim.data.resize(4);

    while( ros::ok() ) {
        
        for( int i=0; i<4; i++ ) {
            f_stim.data[i] = _gamma(i);
        }

        for( int i=0; i<12; i++ ) {
            x_stim.data[i] = _x_tilde[i];
        }

        for( int i=0; i<_motor_num; i++ ) {
            veloc.data[i] = _u_k(i);
        }

        for( int i=0; i<6; i++ ) {
            ykk.data[i] = _y_kk(i);
        }

        _motori.publish( veloc );
        _yy.publish( ykk );
        _x_estim.publish( x_stim );
        _fault_estim.publish( f_stim );
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

void TSKF::tskf_matrix_generation( Eigen::MatrixXd allocation_M ) {
    Eigen::Matrix<double,12,12> I_12;           //matrice identità
    Eigen::Matrix<double,12,12> Ac;
    I_12.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    Ac = _A_k; 
    _A_k = I_12 + Ac*_Ts; //matrice sistema linearizzato discretizzata

    Eigen::Matrix<double,12,4> B;
    B = _B_k;    
    B(5,0) = _motor_force_k/_mass; B(5,1) = _motor_force_k/_mass; B(5,2) = _motor_force_k/_mass; B(5,1) = _motor_force_k/_mass;
    
    B(7,0) = -allocation_M(1,0)/_inertia(1,1); B(7,1) = -allocation_M(1,1)/_inertia(1,1);  
    B(7,2) = -allocation_M(1,2)/_inertia(1,1); B(7,3) = -allocation_M(1,3)/_inertia(1,1);
    
    B(9,0) = allocation_M(0,0)/_inertia(0,0); B(9,1) = allocation_M(0,1)/_inertia(0,0); 
    B(9,2) = allocation_M(0,2)/_inertia(0,0); B(9,3) = allocation_M(0,3)/_inertia(0,0);
    
    B(11,0) = -allocation_M(2,0)/_inertia(2,2); B(11,1) = -allocation_M(2,1)/_inertia(2,2);
    B(11,2) = -allocation_M(2,2)/_inertia(2,2); B(11,3) = -allocation_M(2,3)/_inertia(2,2);

    _B_k = B*_Ts;

    
}

void TSKF::estimation() {
    
    ros::Rate r( 100 );


    //Variabili del tskf usate solo in estimation
    Eigen::Matrix<double,12,4> W_k =Eigen::Matrix<double,12,4>::Zero();             //coupling
    Eigen::Matrix<double,4,4> P_gamma_kk_k = Eigen::Matrix<double,4,4>::Zero();     //predizione P_gamma
    Eigen::Vector4d gamma_pred = Eigen::Vector4d::Zero();                 //predizione gamma
    Eigen::Matrix<double,12,4> V_kk_k = Eigen::Matrix<double,12,4>::Zero();          //coupling
    Eigen::Matrix<double,6,4> H_kk_k = Eigen::Matrix<double,6,4>::Zero();           //coupling
    Eigen::Matrix<double,12,12> P_x_kk_k = Eigen::Matrix<double,12,12>::Zero();       //predizione P_x
    Eigen::Matrix<double,12,1> x_kk_k = Eigen::Matrix<double,12,1>::Zero();          //predizione x
    Eigen::Matrix<double,12,6> Kx_kk = Eigen::Matrix<double,12,6>::Zero();           //guadagno kf per lo stato
    Eigen::Matrix<double,6,6> S_kk = Eigen::Matrix<double,6,6>::Zero();             //residuo
    Eigen::Matrix<double,4,6> K_gamma_kk = Eigen::Matrix<double,4,6>::Zero();       //guadagno kf per gamma
    Eigen::Matrix<double,4,4> U = Eigen::Matrix<double,4,4>::Zero();                //input di controllo (matrice)
    Eigen::Matrix<double,4,4> I_4;              //matrice identità
    Eigen::Matrix<double,12,12> I_12;           //matrice identità

    I_4.diagonal() << 1.0, 1.0, 1.0, 1.0;
    I_12.diagonal() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     1.0, 1.0, 1.0, 1.0, 1.0, 1.0;


    Eigen::MatrixXd allocation_M;
    if(!generate_allocation_matrix( allocation_M, _motor_num, _rotor_angles, _arm_length, _motor_force_k, _motor_moment_k, _motor_rotation_direction ) ) {     
        cout << "Wrong allocation matrix" << endl;
        exit(0);
    }

    tskf_matrix_generation(allocation_M);
    // for ( int i=0; i<12; i++) {
    //     for( int j=0; j<4; j++) {
    //         cout<<_B_k(i,j)<<" , ";
    //     }
    //     cout<<endl;
    // }
        // U << _u_k[0], 0.0, 0.0, 0.0,
        // 0.0, _u_k[1], 0.0, 0.0,
        // 0.0, 0.0, _u_k[2], 0.0,
        // 0.0, 0.0, 0.0, _u_k[3];

        U << 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0;

        W_k = _A_k*_V_kk_kk - _B_k*U; //coupling equation (26)
        // // //W_k =  - _B_k*U;
        P_gamma_kk_k = _P_gamma_kk_kk + _Qgamma; //optimal bias estimator (15)
        gamma_pred = _gamma; //OBE (14)
        V_kk_k = W_k*_P_gamma_kk_kk*P_gamma_kk_k.inverse(); //coupling (27)
        H_kk_k = _C_k*V_kk_k; //coupling (28)
        P_x_kk_k = _A_k*_P_x_kk_kk*_A_k.transpose()+ _Qx  + W_k*_P_gamma_kk_kk*W_k.transpose()
        - V_kk_k*P_gamma_kk_k*V_kk_k.transpose(); //BFE (20)
        x_kk_k = _A_k*_x_tilde + W_k*_gamma - V_kk_k*_gamma + _vec; //BFE (19)
        _res = _y_kk - _C_k*x_kk_k; //calcolo residui (24)
        Kx_kk = P_x_kk_k*_C_k.transpose()*(_C_k*P_x_kk_k*_C_k.transpose() + _R).inverse(); //BFE (22)
        S_kk = _C_k*P_x_kk_k*_C_k.transpose() + _R; //calcolo residui (25)
        K_gamma_kk = P_gamma_kk_k*(H_kk_k.transpose())*(H_kk_k*P_gamma_kk_k*H_kk_k.transpose() + S_kk).inverse();

                //output stimatore
        _V_kk_kk = V_kk_k - Kx_kk*H_kk_k; //coupling (29)
        _P_gamma_kk_kk = (I_4 - K_gamma_kk*H_kk_k)*P_gamma_kk_k; //OBE (18)
        _P_x_kk_kk = (I_12 - Kx_kk*_C_k)*P_x_kk_k; //BFE (23)
        _x_tilde = x_kk_k + Kx_kk*( _y_kk - _C_k*x_kk_k ); //BFE (21)
        _gamma = gamma_pred + K_gamma_kk*( _res - H_kk_k*_gamma); //OBE (16)

        //         for ( int i=0; i<12; i++) {
        //     for( int j=0; j<4; j++) {
        //         cout<< _V_kk_kk(i,j)<<" , ";
        //     }
        //     cout<<endl;
        // }

        // for ( int i=0; i<6; i++) {
        //     cout<<_res(i)<<endl;
        // }   
        // cout<<"fine"<<endl;


    while (ros::ok() ) {

        U << _u_k[0], 0.0, 0.0, 0.0,
        0.0, _u_k[1], 0.0, 0.0,
        0.0, 0.0, _u_k[2], 0.0,
        0.0, 0.0, 0.0, _u_k[3];
        W_k = _A_k*_V_kk_kk - _B_k*U; //coupling equation (26)
        // // //W_k =  - _B_k*U;
        P_gamma_kk_k = _P_gamma_kk_kk + _Qgamma; //optimal bias estimator (15)
        gamma_pred = _gamma; //OBE (14)
        V_kk_k = W_k*_P_gamma_kk_kk*P_gamma_kk_k.inverse(); //coupling (27)
        H_kk_k = _C_k*V_kk_k; //coupling (28)
        P_x_kk_k = _A_k*_P_x_kk_kk*_A_k.transpose() + _Qx + W_k*_P_gamma_kk_kk*W_k.transpose()
        - V_kk_k*P_gamma_kk_k*V_kk_k.transpose(); //BFE (20)
        x_kk_k = _A_k*_x_tilde + W_k*_gamma - V_kk_k*_gamma + _vec; //BFE (19)
        _res = _y_kk - _C_k*x_kk_k; //calcolo residui (24)
        Kx_kk = P_x_kk_k*_C_k.transpose()/*(_C_k*P_x_kk_k*_C_k.transpose() + _R).inverse()*/; //BFE (22)
        S_kk = _C_k*P_x_kk_k*_C_k.transpose() + _R; //calcolo residui (25)
        K_gamma_kk = P_gamma_kk_k*(H_kk_k.transpose())*(H_kk_k*P_gamma_kk_k*H_kk_k.transpose() + S_kk).inverse(); //OBE (17)

        // output stimatore
        _V_kk_kk = V_kk_k - Kx_kk*H_kk_k; //coupling (29)
        _P_gamma_kk_kk = (I_4 - K_gamma_kk*H_kk_k)*P_gamma_kk_k; //OBE (18)
        _P_x_kk_kk = (I_12 - Kx_kk*_C_k)*P_x_kk_k; //BFE (23)
        _x_tilde = x_kk_k + Kx_kk*( _y_kk - _C_k*x_kk_k ); //BFE (21)
        _gamma = gamma_pred + K_gamma_kk*( _res - H_kk_k*_gamma); //OBE (16)

        r.sleep();

        // for ( int i=0; i<6; i++) {
        //     cout<<_res(i)<<endl;
        // }    
        // for ( int i=0; i<12; i++) {
        //     for( int j=0; j<6; j++) {
        //         cout<< Kx_kk(i,j)<<" , ";
        //     }
        //     cout<<endl;
        // }
        // cout<<"fine"<<endl;

    }


}

void TSKF::run() {

    boost::thread estimation_t( &TSKF::estimation, this );
    boost::thread publishr_test_t( &TSKF::publisher_test, this );
    ros::spin();
}

int main( int argc, char** argv ) {

    ros::init(argc, argv, "tskf_estimator");
    TSKF e;
    e.run();

    return 0;
}