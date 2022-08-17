#include <chrono>
#include <iostream>
#include <string>

#include <opencv2/opencv.hpp>
#include <apriltag/apriltag.h>
#include <apriltag/tag36h11.h>
#include <apriltag/tag36h10.h>
#include <apriltag/tag25h9.h>
#include <apriltag/tag16h5.h>
// #include <apriltag/tagCircle21h7.h>
// #include <apriltag/tagCircle49h12.h>
// #include <apriltag/tagCustom48h12.h>
// #include <apriltag/tagStandard41h12.h>
// #include <apriltag/tagStandard52h13.h>


class TagSolver
{
private:
	cv::Mat _intrinsic;
	cv::Mat _dist;
	std::vector<cv::Point3d> _p3d;
	std::vector<cv::Point2d> _p2d;

public:
	TagSolver(double len, double fx, double fy, double cx, double cy) :
		_p3d(4), _p2d(4)
	{
		_intrinsic = cv::Mat(3, 3, CV_64FC1, cv::Scalar(0.0));
		_intrinsic.ptr<double>(0)[0] = fx;
		_intrinsic.ptr<double>(1)[1] = fy;
		_intrinsic.ptr<double>(0)[2] = cx;
		_intrinsic.ptr<double>(1)[2] = cy;
		_intrinsic.ptr<double>(2)[2] = 1.0;

		_dist = cv::Mat(1, 5, CV_64FC1, cv::Scalar(0.0));

		len *= 0.5;
		_p3d[0] = cv::Point3d(-len, len, 0);
		_p3d[1] = cv::Point3d(len, len, 0);
		_p3d[2] = cv::Point3d(len, -len, 0);
		_p3d[3] = cv::Point3d(-len, -len, 0);
	}

	void solve(apriltag::apriltag_detection_t* det, cv::Vec3d& t, cv::Mat1d& rot)
	{
		_p2d[0] = cv::Point2d(det->p[0][0], det->p[0][1]);
		_p2d[1] = cv::Point2d(det->p[1][0], det->p[1][1]);
		_p2d[2] = cv::Point2d(det->p[2][0], det->p[2][1]);
		_p2d[3] = cv::Point2d(det->p[3][0], det->p[3][1]);
		cv::Vec3d r;
		cv::solvePnP(_p3d, _p2d, _intrinsic, _dist, r, t, false, cv::SOLVEPNP_IPPE_SQUARE);
		cv::Rodrigues(r, rot);
	}
};


int main(int argc, char* argv[])
{
	cv::String keys =
		"{h help     |          | print this message}"
		"{d device   | 0        | camera device number}"
		"{iw width   | 0        | width}"
		"{ih height  | 0        | height}"
		"{f family   | tag36h11 | tag family to use}"
		"{x decimate | 2.0      | decimate input image by this factor}"
		"{b blur     | 0.0      | apply low-pass blur to input}"
		"{r refine   | 1        | spend more time trying to align edges of tags}";
	cv::CommandLineParser parser(argc, argv, keys);
	if (parser.has("help"))
	{
		parser.printMessage();
		return 0;
	}

	std::cout << "OpenCV: " << CV_VERSION << std::endl;

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
	// else if (arg_family == "tagCircle21h7")
	// 	tf = apriltag::tagCircle21h7_create();
	// else if (arg_family == "tagCircle49h12")
	// 	tf = apriltag::tagCircle49h12_create();
	// else if (arg_family == "tagStandard41h12")
	// 	tf = apriltag::tagStandard41h12_create();
	// else if (arg_family == "tagStandard52h13")
	// 	tf = apriltag::tagStandard52h13_create();
	// else if (arg_family == "tagCustom48h12")
	// 	tf = apriltag::tagCustom48h12_create();
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
	td->min_white_black_diff = 40;
	td->min_tag_area = 36;
	td->max_line_fit_mse = 1.0;

	const int arg_device = parser.get<int>("device");
	cv::VideoCapture cap(arg_device);
	if (!cap.isOpened())
	{
		std::cout << "Couldn't open video capture device (" << arg_device << ")." << std::endl;
		return -1;
	}
	const int arg_width = parser.get<int>("width");
	if (arg_width > 0)
		cap.set(cv::CAP_PROP_FRAME_WIDTH, arg_width);
	const int arg_height = parser.get<int>("height");
	if (arg_height > 0)
		cap.set(cv::CAP_PROP_FRAME_HEIGHT, arg_height);

	cv::Mat frame, gray;
	while (true)
	{
		cap >> frame;
		cv::cvtColor(frame, gray, cv::COLOR_BGR2GRAY);

		apriltag::image_u8_t im;
		im.width = gray.cols;
		im.height = gray.rows;
		im.buf = gray.data;

		auto beg = std::chrono::steady_clock::now();
		apriltag::zarray_t* detections = apriltag::apriltag_detector_detect(td, &im);
		auto end = std::chrono::steady_clock::now();
		double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(end - beg).count() * 1e-9;
		std::cout << std::endl << "Time detect: " << dt << " s." << std::endl;

		// Draw.
		const cv::Scalar color_top(0, 255, 0);
		const cv::Scalar color(255, 0, 0);
		const int line_w = 2;
		const int text_w = 1;
		for (int i = 0; i < zarray_size(detections); ++i)
		{
			apriltag::apriltag_detection_t* det;
			apriltag::zarray_get(detections, i, &det);
			if (det->hamming > 0)
				continue;

			cv::Point2d pt1(det->p[0][0], det->p[0][1]);
			cv::Point2d pt2(det->p[1][0], det->p[1][1]);
			cv::Point2d pt3(det->p[2][0], det->p[2][1]);
			cv::Point2d pt4(det->p[3][0], det->p[3][1]);
			cv::line(frame, pt2, pt3, color, line_w);
			cv::line(frame, pt3, pt4, color, line_w);
			cv::line(frame, pt4, pt1, color, line_w);
			cv::line(frame, pt1, pt2, color_top, line_w);
			cv::circle(frame, pt1, line_w + 2, color_top, -1);
			cv::circle(frame, pt2, line_w + 2, color, -1);
			cv::circle(frame, pt3, line_w + 2, color, -1);
			cv::circle(frame, pt4, line_w + 2, color, -1);
			int fontface = cv::FONT_HERSHEY_SCRIPT_SIMPLEX;
			double fontscale = 0.7;
			int baseline;
			std::string text = std::to_string(det->id);
			cv::Point c = (pt1 + pt2 + pt3 + pt4) / 4;
			cv::Size textsize = cv::getTextSize(text, fontface, fontscale, 2, &baseline);
			cv::putText(frame, text, cv::Point(c.x - textsize.width / 2, c.y + textsize.height / 2), fontface, fontscale, cv::Scalar(0, 0, 255), 2);

			{
				cv::Vec3d t;
				cv::Mat1d rot(3, 3);
				TagSolver tag_solver(1.0, im.height, im.height, im.width / 2, im.height / 2);
				auto beg = std::chrono::steady_clock::now();
				tag_solver.solve(det, t, rot);
				auto end = std::chrono::steady_clock::now();
				double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(end - beg).count() * 1e-9;
				std::cout << "pose dt = " << dt << std::endl;
				std::cout << "pose = " << t[0] << " " << t[1] << " " << t[2] << std::endl;
			}
		}
		//
		apriltag_detections_destroy(detections);

		cv::imshow("frame", frame);
		int key = cv::waitKey(1);
		if (key == 27)
			break;
		if (key == 32)
			cv::imwrite("frame.png", frame);
	}

	apriltag::apriltag_detector_destroy(td);
	apriltag::tag_destroy(tf);

	return 0;
}