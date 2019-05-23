#include "project_error.h"
#include "pose_local_parameterization.h"
ProjectError::ProjectError(const Eigen::Vector3d& uv_C0):
        C0uv(uv_C0){
}

bool ProjectError::Evaluate(double const *const *parameters,
                                    double *residuals,
                                    double **jacobians) const {
    return EvaluateWithMinimalJacobians(parameters,
                                        residuals,
                                        jacobians, NULL);
}


bool ProjectError::EvaluateWithMinimalJacobians(double const *const *parameters,
                                                        double *residuals,
                                                        double **jacobians,
                                                        double **jacobiansMinimal) const {

    // T_WC
    Eigen::Vector3d t_WC(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Q_WC(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    // Wp
    Eigen::Vector3d Wp(parameters[1][0], parameters[1][1], parameters[1][2]);

    Eigen::Matrix3d R_WC = Q_WC.toRotationMatrix();
    Eigen::Vector3d Cp = R_WC.transpose()*(Wp - t_WC);//将点p从world 转到 camera
    Eigen::Matrix<double, 2, 1> error;

    double inv_z = 1/Cp(2);
    Eigen::Vector2d hat_C0uv(Cp(0)*inv_z, Cp(1)*inv_z);//除法变乘法;得到相机归一化坐标;hat表示观测值???

    Eigen::Matrix<double,2,3> H;
    H << 1, 0, -Cp(0)*inv_z,
            0, 1, -Cp(1)*inv_z;
    H *= inv_z;

    error = hat_C0uv - C0uv.head<2>();//观测值 - 测量值
    squareRootInformation_.setIdentity();
    //squareRootInformation_ = weightScalar_* squareRootInformation_; //Weighted

    // weight it
    Eigen::Map<Eigen::Matrix<double, 2, 1> > weighted_error(residuals);
    weighted_error = squareRootInformation_ * error;



    // calculate jacobians
    if(jacobians != NULL){

        if(jacobians[0] != NULL){
            Eigen::Matrix<double,2,6,Eigen::RowMajor> jacobian0_min;
            Eigen::Map<Eigen::Matrix<double,2,7,Eigen::RowMajor>> jacobian0(jacobians[0]);//引用

            Eigen::Matrix<double, 3, 6> tmp;
            tmp.setIdentity();
            tmp.topLeftCorner(3,3) = - R_WC.transpose();
            tmp.topRightCorner(3,3) =  - Utility::skewSymmetric(R_WC.transpose()*(Wp - t_WC));

            jacobian0_min  =  H*tmp; //?????


            jacobian0 << squareRootInformation_*jacobian0_min, Eigen::Matrix<double, 2,1>::Zero();



            if(jacobiansMinimal != NULL && jacobiansMinimal[0] != NULL){
                Eigen::Map<Eigen::Matrix<double,2,6,Eigen::RowMajor>> map_jacobian0_min(jacobiansMinimal[0]);
                map_jacobian0_min = squareRootInformation_*jacobian0_min;
            }
        }

        if(jacobians[1] != NULL){
            Eigen::Matrix<double,2,3,Eigen::RowMajor> jacobian1_min;
            Eigen::Map<Eigen::Matrix<double,2,3,Eigen::RowMajor>> jacobian1(jacobians[1]);

            Eigen::Matrix<double, 3, 3> tmp;

            tmp = R_WC.transpose();
            jacobian1_min = H*tmp;

            jacobian1 = squareRootInformation_*jacobian1_min;

            if(jacobiansMinimal != NULL && jacobiansMinimal[1] != NULL){
                Eigen::Map<Eigen::Matrix<double,2,3,Eigen::RowMajor>> map_jacobian1_min(jacobiansMinimal[1]);
                map_jacobian1_min = squareRootInformation_*jacobian1_min;
            }
        }

    }
    return true;
}