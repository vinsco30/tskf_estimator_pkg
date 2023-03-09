#include "ros/ros.h"
#include "rosgraph_msgs/Clock.h"
#include "boost/thread.hpp"
#include <Eigen/Eigen>
#include "std_msgs/Float32MultiArray.h"
#include "std_msgs/Bool.h"
#include <Eigen/Geometry>
#include <fstream>
#include <iostream>
//#include "utils.h"

using namespace std;
using namespace Eigen;

class WRITE_RES {

    public:

        WRITE_RES();

        void clock_cb( const rosgraph_msgs::Clock );
        void sys_reset_cb( const std_msgs::Bool );
        void ctrl_act_cb(const std_msgs::Bool );
        void run();
        void write_fun();
    
    private:

        ros::NodeHandle _nh;
        ros::Subscriber _clock_sub;
        ros::Subscriber _sys_res_sub;
        ros::Subscriber _ctrl_sub;

        bool _sys_reset = false;
        bool _ctrl_act = false;
        double _sec = 0;
        double _nsec = 0;
        double test=10;
        int _cnt;


};

WRITE_RES::WRITE_RES() {

    _cnt = 0;
    _clock_sub = _nh.subscribe( "/clock", 0, &WRITE_RES::clock_cb, this );
    _sys_res_sub = _nh.subscribe( "/lee/sys_reset", 0, &WRITE_RES::sys_reset_cb, this );
    _ctrl_sub = _nh.subscribe( "lee/controller_active", 0, &WRITE_RES::ctrl_act_cb, this );

} 

void WRITE_RES::clock_cb( const rosgraph_msgs::Clock t )  {

    _sec = t.clock.toSec();


}

void WRITE_RES::sys_reset_cb( const std_msgs::Bool sr) {

    _sys_reset = sr.data;

}

void WRITE_RES::ctrl_act_cb( const std_msgs::Bool ca) {

    _ctrl_act = ca.data;

}

// void WRITE_RES::write_fun() {
//     std::ofstream csv_file("simulations.csv", std::ios_base::app);

//      csv_file << ros::Time::now().toSec() <<std::endl;
// }

// void WRITE_RES::write_fun() {

//     ros::Rate r(10);

//     // std::fstream data_file;

//     // data_file.open ("simulations.csv");
//     std::ofstream csv_file("simulations.csv", std::ios_base::app);
//     // if( data_file.is_open()) {
//     //     ROS_INFO("File open");
//     // }
//     // else {
//     //     ROS_INFO("Error opening file");
//     //     exit(0);
//     // }

//     data_file << "This is the first cell in the first column.\n";
//     data_file << "a,b,c,\n";
//     data_file << "c,s,v,\n";
//     data_file << "1,2,3.456\n";
//     data_file << "semi;colon";
//         ROS_INFO("File close");
//     data_file.close();
//     r.sleep();



// }

void WRITE_RES::run() {

    //boost::thread write_t( &WRITE_RES::write_fun, this );
    ros::spin();
}

int main( int argc, char** argv ) {

    ros::init(argc, argv, "save_data");
    WRITE_RES s;
    s.run();

    return 0;
}



