# ud-meeg671-final-project
Create joint angle trajectory for a 7 degree-of-freedom robot arm (LBR iiwa 7 R800, KUKA) in MATLAB to make the end effector reach the target.

A 7 degree-of-freedom robot arm is mounted on a base. A wooden flange is mounted on the robot flange. The wooden flange holds a camera. When the camera is not used, the wooden flange is magnetically holding a wooden object. The goal of this project is to program a trajectory for the robot to position the wooden object on top of a rectangular target located in the workspace of the robot arm. The exact location and orientation of the target can be computed using an Aruco marker, that is printed close to the target, and is seen by the camera, when the camera is attachedon the wooden flange. The robot starts from a given confguration q1 and needs to position the wooden object as close as possible to the target, matching its corners A, B, C, D with those of the target. During the motion,the robot should not hit in any way the surface/table the target is put on.

Design the trajectory in joint space so that the robot arm moves from its starting configuration and finishes when the wooden object is on top of the target. In order to locate the target position, the camera takes a picture of the target and the Aruco marker. The software returns the position and orientation of the Aruco marker reference system with respect to the camera reference system. When creating the robot trajectory make sure the angular position and angular velocity limits of all the robot joints don't exceed. The entire robot motion should not last for more than 30 seconds.
