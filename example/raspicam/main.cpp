#include <chrono>
#include <iostream>
#include <string>
#include <memory>

// https://github.com/cedricve/raspicam
#include <raspicam/raspicam.h>
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


int main(int argc, char* argv[])
{
	cv::String keys =
		"{h help     |          | print this message}"
		"{iw width   | 1280     | width}"
		"{ih height  | 1024     | height}"
		"{fps        | 30       | fps}"
		"{iso        | 800      | iso}"
		"{f family   | tag36h11 | tag family to use}"
		"{x decimate | 2.0      | decimate input image by this factor}"
		"{b blur     | 0.0      | apply low-pass blur to input}"
		"{r refine   | 1        | spend more time trying to align edges of tags}"
		"{s scale    | 1.0      | image scale}";
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

	const int arg_width = parser.get<int>("width");
	const int arg_height = parser.get<int>("height");
	const int arg_fps = parser.get<int>("fps");
	const int arg_iso = parser.get<int>("iso");
	double arg_scale = parser.get<double>("scale");
	if (arg_scale < 0.1)
		arg_scale = 0.1;

	raspicam::RaspiCam cam;
	cam.setFormat(raspicam::RASPICAM_FORMAT_GRAY);
	cam.setCaptureSize(arg_width, arg_height);
	cam.setFrameRate(arg_fps);
	cam.setISO(arg_iso);

	int width = cam.getWidth();
	int height = cam.getHeight();
	if (!cam.open())
	{
		std::cout << "camera not opended." << std::endl;
		return 0;
	}
	std::cout << "raspicam: " << width << "x" << height << "@" << cam.getFrameRate() << std::endl;

	cv::Mat frame, gray_sized;;
	while (true)
	{
		auto beg = std::chrono::steady_clock::now();
		if (!cam.grab())
			continue;

		apriltag::image_u8_t im;
		im.width = width;
		im.height = height;
		im.buf = cam.getImageBufferData();

		apriltag::zarray_t* detections = apriltag::apriltag_detector_detect(td, &im);
		auto end = std::chrono::steady_clock::now();
		double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(end - beg).count() * 1e-9;
		std::cout << "Time detect: " << dt << " s." << std::endl;

		// Draw.
		cv::Mat gray(height, width, CV_8UC1, cam.getImageBufferData());
		if (arg_scale == 1.0)
			gray_sized = gray;
		else
			cv::resize(gray, gray_sized, cv::Size(), arg_scale, arg_scale);
		cv::cvtColor(gray_sized, frame, cv::COLOR_GRAY2BGR);
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
			pt1 *= arg_scale;
			pt2 *= arg_scale;
			pt3 *= arg_scale;
			pt4 *= arg_scale;
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

	cam.release();

	return 0;
}