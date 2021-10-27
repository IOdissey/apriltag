#include <chrono>
#include <iostream>
#include <string>

#include <opencv2/opencv.hpp>
#include <ros/ros.h>
#include <sensor_msgs/CompressedImage.h>
#include <apriltag/apriltag.h>
#include <apriltag/tag36h11.h>
#include <apriltag/tag36h10.h>
#include <apriltag/tag25h9.h>
#include <apriltag/tag16h5.h>


cv::Mat frame;
bool is_frame = false;
void callback_cam(const sensor_msgs::CompressedImageConstPtr& msg)
{
	int size = msg->data.size();
	const uchar* data = msg->data.data();
	cv::Mat raw(1, size, CV_8UC1, (void*)data);
	frame = cv::imdecode(raw, 1);
	is_frame = true;
}

int main(int argc, char* argv[])
{
	cv::String keys =
		"{h help        |          | print this message}"
		"{t topic       | 0        | ros topic}"
		"{f family      | tag16h5  | tag family to use}"
		"{x decimate    | 2.0      | decimate input image by this factor}"
		"{b blur        | 0.0      | apply low-pass blur to input}"
		"{r refine      | 1        | spend more time trying to align edges of tags}"
		"{bd black_diff | 50       | }";
	cv::CommandLineParser parser(argc, argv, keys);
	if (parser.has("help"))
	{
		parser.printMessage();
		return 0;
	}

	apriltag::apriltag_family_t* tf = nullptr;
	const std::string arg_family = parser.get<std::string>("family");
	if (arg_family == "tag36h11")
		tf = apriltag::tag36h11_create();
	else if (arg_family == "tag36h10")
		tf = apriltag::tag36h10_create();
	else if (arg_family == "tag25h9")
		tf = apriltag::tag25h9_create();
	else if (arg_family == "tag16h5")
		tf = apriltag::tag16h5_create();
	else
	{
		std::cout << "Unrecognized tag family name (" << arg_family << ")." << std::endl;
		return -1;
	}

	apriltag::apriltag_detector_t* td = apriltag::apriltag_detector_create();
	apriltag::apriltag_detector_add_family(td, tf);
	td->quad_decimate = parser.get<float>("decimate");
	td->quad_sigma = parser.get<float>("blur");
	td->refine_edges = parser.get<int>("refine");
	td->qtp.min_white_black_diff = parser.get<int>("black_diff");

	// Запуск ROS.
	ros::init(argc, argv, "apriltag_ros");
	ros::Time::init();
	if (!ros::ok())
	{
		std::cout << "ROS init error." << std::endl;
		return -1;
	}
	if ((!ros::master::check()))
	{
		std::cout << "ROS master not started." << std::endl;
		return -1;
	}
	ros::NodeHandle node("~");
	ros::Subscriber sub_cam = node.subscribe(parser.get<std::string>("topic"), 1, &callback_cam);
	ros::Rate rate(100);
	std::cout << "ROS topic: " << parser.get<std::string>("topic") << std::endl;

	cv::Mat gray;
	while (true)
	{
		if (!ros::ok())
			break;
		ros::spinOnce();
		if (is_frame)
		{
			is_frame = false;
			cv::cvtColor(frame, gray, CV_BGR2GRAY);

			apriltag::image_u8_t im;
			im.width = gray.cols;
			im.height = gray.rows;
			im.buf = gray.data;

			auto beg = std::chrono::steady_clock::now();
			apriltag::zarray_t* detections = apriltag::apriltag_detector_detect(td, &im);
			auto end = std::chrono::steady_clock::now();
			double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(end - beg).count() * 1e-9;
			std::cout << "Time detect: " << dt << " s." << std::endl;

			// Draw.
			const cv::Scalar color(0, 255, 0);
			const int line_w = 3;
			const int text_w = 3;
			for (int i = 0; i < zarray_size(detections); ++i)
			{
				apriltag::apriltag_detection_t* det;
				apriltag::zarray_get(detections, i, &det);
				cv::Point2d pt1(det->p[0][0], det->p[0][1]);
				cv::Point2d pt2(det->p[1][0], det->p[1][1]);
				cv::Point2d pt3(det->p[2][0], det->p[2][1]);
				cv::Point2d pt4(det->p[3][0], det->p[3][1]);
				cv::line(frame, pt1, pt2, color, line_w);
				cv::line(frame, pt2, pt3, color, line_w);
				cv::line(frame, pt3, pt4, color, line_w);
				cv::line(frame, pt4, pt1, color, line_w);
				int fontface = cv::FONT_HERSHEY_SCRIPT_SIMPLEX;
				double fontscale = 0.7;
				int baseline;
				std::string text = std::to_string(det->id);
				cv::Point c = (pt1 + pt2 + pt3 + pt4) / 4;
				cv::Size textsize = cv::getTextSize(text, fontface, fontscale, 2, &baseline);
				cv::putText(frame, text, cv::Point(c.x - textsize.width / 2, c.y + textsize.height / 2), fontface, fontscale, cv::Scalar(0, 0, 255), 2);
			}

			cv::imshow("frame", frame);
			int key = cv::waitKey(1);
			if (key == 27)
				break;
			if (key == 32)
				cv::imwrite("frame.png", frame);
		}
		rate.sleep();
	}

	apriltag::apriltag_detector_destroy(td);
	apriltag::tag_destroy(tf);

	return 0;
}