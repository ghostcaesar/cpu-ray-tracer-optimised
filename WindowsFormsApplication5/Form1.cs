using System;
using System.Drawing;
using System.Windows.Forms;
using System.Threading;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;

namespace WindowsFormsApplication5
{



    public partial class Form1 : Form
    {
        bool th1 = false;
        bool th2 = false;
        bool th3 = false;
        bool th4 = false;

        void paint()
        {
            output.Image = Picture.bmpPic;
        }

        public Form1()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {

        }

        private void Form1_Load(object sender, EventArgs e)
        {

            //thread 1
            Worker worker1 = new Worker();
            Thread thd1 = new Thread(worker1.DoWork1);

            //thread 2
            Thread thd2 = new Thread(worker1.DoWork2);

            //thread 1
            Thread thd3 = new Thread(worker1.DoWork3);

            //thread 4
            Thread thd4 = new Thread(worker1.DoWork4);

            thd1.Start();
            thd2.Start();
            thd3.Start();
            thd4.Start();


            while (thd1.IsAlive || thd2.IsAlive || thd3.IsAlive || thd4.IsAlive)
            {
                Thread.Sleep(20);
            }
            Thread.Sleep(100);
            //Picture.setBmp();
            paint();
                
            Picture.bmpPic.Save("myBitmap.bmp");

        }


    }

    class Worker
    {
        public void DoWork1()
        {
            Picture pic1 = new Picture();
            EIn objIn = new EIn();
            objIn.start = 00;
            objIn.end = 1001000;
            pic1.render(objIn);
        }

        public void DoWork2()
        {
            Picture pic1 = new Picture();
            EIn objIn = new EIn();
            objIn.start = 1001001;
            objIn.end = 2002000;
            pic1.render(objIn);
        }

        public void DoWork3()
        {
            Picture pic1 = new Picture();
            EIn objIn = new EIn();
            objIn.start = 2002001;
            objIn.end = 3003000;
            pic1.render(objIn);
        }

        public void DoWork4()
        {
            Picture pic1 = new Picture();
            EIn objIn = new EIn();
            objIn.start = 3003001;
            objIn.end = 4004000;
            pic1.render(objIn);
        }
    }

    public class EIn
    {
        public int start;
        public int end;
    }

    public class Sphere
    {


        public float[] pos = { 0, 0, 0 };

        public float radius = 0;

        //reflection, refraction
        public float[] attributes = { 0.25f, 0 };

        public float[] colour = { 0, 0, 0 };

        public Sphere(float x, float y, float z, float r, float[] rgb)
        {
            pos[0] = x;
            pos[1] = y;
            pos[2] = z;
            radius = r;
            colour = rgb;
        }
        public object clone()
        {
            return this.MemberwiseClone();
        }
    }

    public class Plane
    {


        public double[] normal = { 0, 0, 0 };

        public double dist = 0;

        //reflection, refraction
        public float[] attributes = { 0.7f, 0 };

        public float[] colour = { 0, 0, 0 };

        public Plane(double normal_x, double normal_y, double normal_z, double dist, float[] rgb)
        {
            normal[0] = normal_x;
            normal[1] = normal_y;
            normal[2] = normal_z;
            this.dist = dist;
            colour = rgb;
        }

        public object clone()
        {
            return this.MemberwiseClone();
        }
    }

    class PointLight
    {
        public float[] pos = { -170, 440, 500 };

        public float intensity = 100;
    }

    class Scene
    {
        public int reflectDepth = 4;
        public PointLight light = new PointLight();

        // 1/255 of RGB val
        private static float[] sph1_colour = { 0.7f, 0, 0 };
        private static float[] sph2_colour = { 1, 1, 1 };

        public Sphere[] spheres = new Sphere[]
        {
            new Sphere(0, 420, 70, 40,sph1_colour),
            new Sphere(-40, 420, -70,40,sph2_colour),
        };

        // 1/255 of RGB val
        private static float[] p1_colour = { 0.2f, 0.2f, 0.5F };
        private static float[] p2_colour = { 0.8f, 0.8f, 0.8f };
        private static float[] p3_colour = { 0.2f, 0.2f, 0.4f };

        public Plane[] planes = new Plane[]
        {
            new Plane(1,0,0, 125,p1_colour),
            new Plane(0,0,1,-110,p2_colour),
            new Plane(0,1,0,500,p3_colour),
        };

    }

    class Picture
    {
        Scene scene = new Scene();

        
        private static byte[,,] pixelb = new byte[2001, 2001, 4];

        private static IntPtr pointer;

        public static Bitmap bmpPic;

        static PixelFormat pxf1 = PixelFormat.Format32bppRgb;

        public Picture()
        {
            GCHandle handle = GCHandle.Alloc(pixelb, GCHandleType.Pinned);
            try
            {
                pointer = handle.AddrOfPinnedObject();
            }
            finally
            {
                if (handle.IsAllocated)
                {
                    handle.Free();
                }
            }

            bmpPic = new Bitmap(2001, 2001, 8004,pxf1, pointer);

    }





        float[] cameraPos = { 0, 0, 0 };

        //position of screen centre
        float[] screenPos = { 0, 1, 0 };
        //position of unit vector for screen direction
        float[] screenDir = { 0, 1, 0 };


        private double[] sphereIntersect(int x_count, int y_count, Sphere s_in)
        {
            double[] ray_dir = new double[3];
            ray_dir[0] = x_count; //temp
            ray_dir[1] = 1000; //temp
            ray_dir[2] = y_count; //temp

            double unit_mult = 1 / Math.Sqrt(ray_dir[0] * ray_dir[0] + ray_dir[1] * ray_dir[1] + ray_dir[2] * ray_dir[2]);

            ray_dir[0] = ray_dir[0] * unit_mult;
            ray_dir[1] = ray_dir[1] * unit_mult;
            ray_dir[2] = ray_dir[2] * unit_mult;

            //surface point co-ord, 1 t_length, 3 ray direction
            double[] ret = { 0, 0, 0, 0, ray_dir[0], ray_dir[1], ray_dir[2] };

            double a = ray_dir[0] * ray_dir[0] + ray_dir[1] * ray_dir[1] + ray_dir[2] * ray_dir[2];
            double b = -2 * (ray_dir[0] * s_in.pos[0] + ray_dir[1] * s_in.pos[1] + ray_dir[2] * s_in.pos[2]);
            double c = s_in.pos[0] * s_in.pos[0] + s_in.pos[1] * s_in.pos[1] + s_in.pos[2] * s_in.pos[2] - s_in.radius * s_in.radius;

            //b*b -4ac
            double d = b * b - 4 * a * c;
            if (d < 0)
                return ret;

            double t1 = -0.5 * (b - Math.Sqrt(d)) / a;
            double t2 = -0.5 * (b + Math.Sqrt(d)) / a;

            if (isPositive(t1))
            {
                if (isPositive(t2))
                {
                    if (t1 < t2)
                    {
                        ret[0] = t1 * ray_dir[0];
                        ret[1] = t1 * ray_dir[1];
                        ret[2] = t1 * ray_dir[2];
                        ret[3] = t1;
                    }
                    else
                    {
                        ret[0] = t2 * ray_dir[0];
                        ret[1] = t2 * ray_dir[1];
                        ret[2] = t2 * ray_dir[2];
                        ret[3] = t2;
                    }
                }
                else
                {
                    ret[0] = t1 * ray_dir[0];
                    ret[1] = t1 * ray_dir[1];
                    ret[2] = t1 * ray_dir[2];
                    ret[3] = t1;
                }
            }
            else if (isPositive(t2))
            {
                ret[0] = t2 * ray_dir[0];
                ret[1] = t2 * ray_dir[1];
                ret[2] = t2 * ray_dir[2];
                ret[3] = t2;
            }


            return ret;
        }



        private double[] planeIntersect(int x_count, int y_count, Plane p_in)
        {

            double[] ray_dir = new double[3];
            ray_dir[0] = x_count; //temp
            ray_dir[1] = 1000; //temp
            ray_dir[2] = y_count; //temp

            double unit_mult = 1 / Math.Sqrt(ray_dir[0] * ray_dir[0] + ray_dir[1] * ray_dir[1] + ray_dir[2] * ray_dir[2]);

            ray_dir[0] = ray_dir[0] * unit_mult;
            ray_dir[1] = ray_dir[1] * unit_mult;
            ray_dir[2] = ray_dir[2] * unit_mult;

            //3 surface point, 1 t_length, 3 ray direction
            double[] ret = { 0, 0, 0, 0, ray_dir[0], ray_dir[1], ray_dir[2] };

            double t1 = (p_in.dist - 0) / (ray_dir[0] * p_in.normal[0] + ray_dir[1] * p_in.normal[1] + ray_dir[2] * p_in.normal[2]);

            if (isPositive(t1))
            {
                ret[0] = t1 * ray_dir[0];
                ret[1] = t1 * ray_dir[1];
                ret[2] = t1 * ray_dir[2];
                ret[3] = t1;
            }

            return ret;
        }

        private double[] FindReflectedRayPlane(double[] sur_in, Plane pl_in)
        {
            double[] ret = new double[3];

            double[] normal = new double[3];

            normal[0] = pl_in.normal[0];
            normal[1] = pl_in.normal[1];
            normal[2] = pl_in.normal[2];


            //comment out for performance, input normal should already be unit size
            //double normal_unit_multiplier = 1 / Math.Sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

            ////unit normal
            //normal[0] = normal_unit_multiplier * normal[0];
            //normal[1] = normal_unit_multiplier * normal[1];
            //normal[2] = normal_unit_multiplier * normal[2];


            double n_dot_d = normal[0] * sur_in[4] + normal[1] * sur_in[5] + normal[2] * sur_in[6];
            ret[0] = sur_in[4] - 2 * n_dot_d * normal[0];
            ret[1] = sur_in[5] - 2 * n_dot_d * normal[1];
            ret[2] = sur_in[6] - 2 * n_dot_d * normal[2];

            double ret_unit_multiplier = 1 / Math.Sqrt(ret[0] * ret[0] + ret[1] * ret[1] + ret[2] * ret[2]);
            ret[0] = ret[0] * ret_unit_multiplier;
            ret[1] = ret[1] * ret_unit_multiplier;
            ret[2] = ret[2] * ret_unit_multiplier;


            return ret;
        }

        private double[] FindReflectedRaySph(double[] sur_in, Sphere sph_in)
        {
            double[] ret = new double[3];

            double[] normal = new double[3];

            normal[0] = sur_in[0] - sph_in.pos[0];
            normal[1] = sur_in[1] - sph_in.pos[1];
            normal[2] = sur_in[2] - sph_in.pos[2];

            double normal_unit_multiplier = 1 / Math.Sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

            //unit normal
            normal[0] = normal_unit_multiplier * normal[0];
            normal[1] = normal_unit_multiplier * normal[1];
            normal[2] = normal_unit_multiplier * normal[2];


            double n_dot_d = normal[0] * sur_in[4] + normal[1] * sur_in[5] + normal[2] * sur_in[6];
            ret[0] = sur_in[4] - 2 * n_dot_d * normal[0];
            ret[1] = sur_in[5] - 2 * n_dot_d * normal[1];
            ret[2] = sur_in[6] - 2 * n_dot_d * normal[2];

            double ret_unit_multiplier = 1 / Math.Sqrt(ret[0] * ret[0] + ret[1] * ret[1] + ret[2] * ret[2]);
            ret[0] = ret[0] * ret_unit_multiplier;
            ret[1] = ret[1] * ret_unit_multiplier;
            ret[2] = ret[2] * ret_unit_multiplier;


            return ret;
        }

        private double[] reflectionRecursiveSph(double[] t_in, Sphere s_in, double[] ray_dir)
        {

            double[] ret = { 0, 0, 0, 0, ray_dir[0], ray_dir[1], ray_dir[2] };

            double a = ray_dir[0] * ray_dir[0] + ray_dir[1] * ray_dir[1] + ray_dir[2] * ray_dir[2];
            double b = -2 * (ray_dir[0] * (s_in.pos[0] - t_in[0]) + ray_dir[1] * (s_in.pos[1] - t_in[1]) + ray_dir[2] * (s_in.pos[2] - t_in[2]));//check sign
            double c = (s_in.pos[0] - t_in[0]) * (s_in.pos[0] - t_in[0]) + (s_in.pos[1] - t_in[1]) * (s_in.pos[1] - t_in[1]) + (s_in.pos[2] - t_in[2]) * (s_in.pos[2] - t_in[2]) - s_in.radius * s_in.radius;

            double d = b * b - 4 * a * c;
            if (d < 0)
                return ret;

            double t1 = -0.5 * (b - Math.Sqrt(d)) / a;
            double t2 = -0.5 * (b + Math.Sqrt(d)) / a;

            if (isPositive(t1) && isPositive(t2))
            {
                if (t1 < t2)
                {
                    ret[0] = t1 * ray_dir[0] + t_in[0];
                    ret[1] = t1 * ray_dir[1] + t_in[1];
                    ret[2] = t1 * ray_dir[2] + t_in[2];
                    ret[3] = t1;
                }
                else
                {
                    ret[0] = t2 * ray_dir[0] + t_in[0];
                    ret[1] = t2 * ray_dir[1] + t_in[1];
                    ret[2] = t2 * ray_dir[2] + t_in[2];
                    ret[3] = t2;
                }
            }
            else if (isPositive(t1))
            {
                ret[0] = t1 * ray_dir[0] + t_in[0];
                ret[1] = t1 * ray_dir[1] + t_in[1];
                ret[2] = t1 * ray_dir[2] + t_in[2];
                ret[3] = t1;
            }
            else if (isPositive(t2))
            {
                ret[0] = t2 * ray_dir[0] + t_in[0];
                ret[1] = t2 * ray_dir[1] + t_in[1];
                ret[2] = t2 * ray_dir[2] + t_in[2];
                ret[3] = t2;
            }


            return ret;
        }

        private double[] reflectionRecursivePlane(double[] t_in, Plane pl_in, double[] ray_dir)
        {

            double[] ret = { 0, 0, 0, 0, ray_dir[0], ray_dir[1], ray_dir[2] };

            double offsetStart = t_in[0] * pl_in.normal[0] + t_in[1] * pl_in.normal[1] + t_in[2] * pl_in.normal[2];
            double t1 = (pl_in.dist - offsetStart - 0) / (ray_dir[0] * pl_in.normal[0] + ray_dir[1] * pl_in.normal[1] + ray_dir[2] * pl_in.normal[2]);

            if (isPositive(t1))
            {
                ret[0] = t1 * ray_dir[0] + t_in[0];
                ret[1] = t1 * ray_dir[1] + t_in[1];
                ret[2] = t1 * ray_dir[2] + t_in[2];
                ret[3] = t1;
            }
            return ret;

        }

        //direct ray from intersect to light source
        //gives intensity if clear
        //else give 0,0,0
        private bool ShadowRaySph(double[] t_in, Sphere s_in)
        {
            double[] ray_dir = new double[3];
            ray_dir[0] = -t_in[0] + scene.light.pos[0]; //temp
            ray_dir[1] = -t_in[1] + scene.light.pos[1]; //temp
            ray_dir[2] = -t_in[2] + scene.light.pos[2]; //temp

            float[] ret = { 0, 0, 0 };

            double a = ray_dir[0] * ray_dir[0] + ray_dir[1] * ray_dir[1] + ray_dir[2] * ray_dir[2];
            //double b = -2 * (ray_dir[0] * (s_in.pos[0]) + ray_dir[1] * (s_in.pos[1]) + ray_dir[2] * (s_in.pos[2] ));//check sign

            double b = -2 * (ray_dir[0] * (s_in.pos[0] - t_in[0]) + ray_dir[1] * (s_in.pos[1] - t_in[1]) + ray_dir[2] * (s_in.pos[2] - t_in[2]));//check sign
            double c = (s_in.pos[0] - t_in[0]) * (s_in.pos[0] - t_in[0]) + (s_in.pos[1] - t_in[1]) * (s_in.pos[1] - t_in[1]) + (s_in.pos[2] - t_in[2]) * (s_in.pos[2] - t_in[2]) - s_in.radius * s_in.radius;

            double d = b * b - 4 * a * c;

            if (d < 0)
                return true;

            double t1 = (-b - Math.Sqrt(d)) / (2 * a);
            double t2 = (-b + Math.Sqrt(d)) / (2 * a);

            if (t1 < -0.00001 || t2 < -0.00001)// take into account NaN
            {
                return true;
                //ret[0] = scene.light.intensity * s_in.attributes[0];
                //ret[1] = scene.light.intensity * s_in.attributes[0];
                //ret[2] = scene.light.intensity * s_in.attributes[0];
            }

            return false;
        }

        private bool ShadowRayPlane(double[] t_in, Plane pl_in)
        {

            double[] ray_dir = new double[3];
            ray_dir[0] = -t_in[0] + scene.light.pos[0]; //temp
            ray_dir[1] = -t_in[1] + scene.light.pos[1]; //temp
            ray_dir[2] = -t_in[2] + scene.light.pos[2]; //temp

            float[] ret = { 0, 0, 0 };

            double offsetStart = t_in[0] * pl_in.normal[0] + t_in[1] * pl_in.normal[1] + t_in[2] * pl_in.normal[2];
            double t1 = (pl_in.dist - offsetStart) / (ray_dir[0] * pl_in.normal[0] + ray_dir[1] * pl_in.normal[1] + ray_dir[2] * pl_in.normal[2]);

            if (t1 > 0.00001 && t1 < 0.9999)// take into account NaN
            {
                return false;
                //ret[0] = scene.light.intensity * s_in.attributes[0];
                //ret[1] = scene.light.intensity * s_in.attributes[0];
                //ret[2] = scene.light.intensity * s_in.attributes[0];
            }

            return true;
        }

        private bool isNearZero(double inputNumber)
        {
            if (inputNumber < 0.00001 && inputNumber > -0.00001)
                return true;
            else
                return false;
        }

        private bool isPositive(double inputNumber)
        {
            if (inputNumber > 0.00001)
                return true;

            return false;
        }

        public void render(EIn objIn)
        {

            //screen rotation()
            //screen translation()

            float[] tmp_colour0 = { 0, 0, 0 };
            Sphere sph_tmp = new Sphere(0, 0, 0, -1, tmp_colour0);
            Plane plane_tmp = new Plane(0, 0, 0, -1, tmp_colour0);
            Sphere sph_tmp2 = new Sphere(0, 0, 0, -1, tmp_colour0);
            Plane plane_tmp2 = new Plane(0, 0, 0, -1, tmp_colour0);

            int startx = objIn.start % 2001;
            int starty = objIn.start / 2001;
            int endx = objIn.end % 2001;
            int endy = objIn.end / 2001;

            //camera translation()
            for (int y_dir = starty - 1000; y_dir <= endy - 1000; y_dir++)
            {
                int row_start, row_end;

                if (y_dir == starty)
                    row_start = startx;
                else
                    row_start = 0;

                if (y_dir == endy)
                    row_end = endx;
                else
                    row_end = 2000;


                for (int x_dir = row_start - 1000; x_dir <= row_end - 1000; x_dir++)
                {
                    float intensity = 1;
                    float reflectionValue = 1;

                    bool temp_pixel_val = true;
                    int reflectDepth = 0;
                    float lightingMultiplier;
                    double[] reflectedRay = { 0, 0, 0 };

                    //ambient Lighting
                    float[] pixelColour = { 0, 0, 0 };

                    double[] t_in = { 0, 0, 0, 0, 0, 0, 0 };

                    double[] t_in_tmp = { 0, 0, 0, 0, 0, 0, 0 };

                    foreach (Sphere sph in scene.spheres)
                    {
                        t_in_tmp = sphereIntersect(x_dir, y_dir, sph);
                        if (t_in_tmp[3] != 0 && ((t_in_tmp[3] < t_in[3]) || t_in[3] == 0))
                        {
                            Array.Copy(t_in_tmp, t_in, 7);

                            sph_tmp = (Sphere)sph.clone();
                            //1st is not plane
                            plane_tmp.dist = -1;
                        }
                    }
                    foreach (Plane plane in scene.planes)
                    {
                        t_in_tmp = planeIntersect(x_dir, y_dir, plane);
                        if (t_in_tmp[3] != 0 && ((t_in_tmp[3] < t_in[3]) || t_in[3] == 0))
                        {
                            Array.Copy(t_in_tmp, t_in, 7);

                            plane_tmp = (Plane)plane.clone();
                            //1st is not sphere
                            sph_tmp.radius = -1;
                        }
                    }

                    if (t_in[3] != 0)
                    {
                        intensity = (float)Math.Pow(0.9995, t_in[3] + 1);
                        foreach (Sphere sph in scene.spheres)
                        {
                            temp_pixel_val &= ShadowRaySph(t_in, sph);

                        }
                        foreach (Plane plane in scene.planes)
                        {
                            temp_pixel_val &= ShadowRayPlane(t_in, plane);
                        }
                        //no obstruction between first intersection and light source
                        if (temp_pixel_val)
                        {
                            //if the first object intersected is a sphere
                            if (sph_tmp.radius != -1)
                            {

                                lightingMultiplier = scene.light.intensity * intensity;
                                pixelColour[0] += lightingMultiplier * sph_tmp.colour[0];
                                pixelColour[1] += lightingMultiplier * sph_tmp.colour[1];
                                pixelColour[2] += lightingMultiplier * sph_tmp.colour[2];
                            }
                            //if the first object intersected is a plane
                            else if (plane_tmp.dist != -1)
                            {
                                lightingMultiplier = scene.light.intensity * intensity;
                                pixelColour[0] += lightingMultiplier * plane_tmp.colour[0];
                                pixelColour[1] += lightingMultiplier * plane_tmp.colour[1];
                                pixelColour[2] += lightingMultiplier * plane_tmp.colour[2];
                            }

                        }
                        //obstruction between first intersection and light source
                        else
                        {
                            if (sph_tmp.radius != -1)
                            {
                                lightingMultiplier = scene.light.intensity * 0.65f * intensity;
                                pixelColour[0] += lightingMultiplier * sph_tmp.colour[0];
                                pixelColour[1] += lightingMultiplier * sph_tmp.colour[1];
                                pixelColour[2] += lightingMultiplier * sph_tmp.colour[2];
                            }
                            else if (plane_tmp.dist != -1)
                            {
                                lightingMultiplier = scene.light.intensity * 0.65f * intensity;
                                pixelColour[0] += lightingMultiplier * plane_tmp.colour[0];
                                pixelColour[1] += lightingMultiplier * plane_tmp.colour[1];
                                pixelColour[2] += lightingMultiplier * plane_tmp.colour[2];
                            }
                        }
                        temp_pixel_val = true;
                    }

                    while (reflectDepth < 3)
                    {
                        double[] t_in_new = { 0, 0, 0, 0, 0, 0, 0 };

                        if (sph_tmp.radius != -1)
                            reflectedRay = FindReflectedRaySph(t_in, sph_tmp);
                        else if (plane_tmp.dist != -1)
                            reflectedRay = FindReflectedRayPlane(t_in, plane_tmp);


                        foreach (Sphere sph in scene.spheres)
                        {

                            t_in_tmp = reflectionRecursiveSph(t_in, sph, reflectedRay);


                            if (t_in_tmp[3] != 0 && (t_in_tmp[3] < t_in_new[3] || t_in_new[3] == 0))
                            {
                                Array.Copy(t_in_tmp, t_in_new, 7);

                                sph_tmp2 = (Sphere)sph.clone();

                                plane_tmp2.dist = -1;
                            }
                        }
                        foreach (Plane plane in scene.planes)
                        {
                            t_in_tmp = reflectionRecursivePlane(t_in, plane, reflectedRay);

                            if (t_in_tmp[3] != 0 && (t_in_tmp[3] < t_in_new[3] || t_in_new[3] == 0))
                            {
                                Array.Copy(t_in_tmp, t_in_new, 7);

                                plane_tmp2 = (Plane)plane.clone();

                                sph_tmp2.radius = -1;
                            }
                        }

                        Array.Copy(t_in_new, t_in, 7);

                        for (int i = 0; i < 7; i++)
                            t_in_new[i] = 0;

                        plane_tmp = (Plane)plane_tmp2.clone();
                        sph_tmp = (Sphere)sph_tmp2.clone();

                        plane_tmp2.dist = -1;
                        sph_tmp2.radius = -1;


                        if (!isNearZero(t_in[3]))
                        {
                            intensity = intensity * (float)Math.Pow(0.9995, t_in[3] + 1);
                            foreach (Sphere sph in scene.spheres)
                            {
                                temp_pixel_val &= ShadowRaySph(t_in, sph);
                            }
                            foreach (Plane plane in scene.planes)
                            {
                                temp_pixel_val &= ShadowRayPlane(t_in, plane);
                            }
                            if (temp_pixel_val)
                            {
                                if (sph_tmp.radius != -1)
                                {
                                    reflectionValue = reflectionValue * sph_tmp.attributes[0];
                                    lightingMultiplier = scene.light.intensity * reflectionValue * intensity;
                                    pixelColour[0] += lightingMultiplier * sph_tmp.colour[0];
                                    pixelColour[1] += lightingMultiplier * sph_tmp.colour[1];
                                    pixelColour[2] += lightingMultiplier * sph_tmp.colour[2];
                                }
                                else if (plane_tmp.dist != -1)
                                {
                                    reflectionValue = reflectionValue * plane_tmp.attributes[0];
                                    lightingMultiplier = scene.light.intensity * reflectionValue * intensity;
                                    pixelColour[0] += lightingMultiplier * plane_tmp.colour[0];
                                    pixelColour[1] += lightingMultiplier * plane_tmp.colour[1];
                                    pixelColour[2] += lightingMultiplier * plane_tmp.colour[2];
                                }
                            }
                            else
                            {
                                if (sph_tmp.radius != -1)
                                {
                                    reflectionValue = reflectionValue * sph_tmp.attributes[0];
                                    lightingMultiplier = scene.light.intensity * reflectionValue * 0.65f * intensity;
                                    pixelColour[0] += lightingMultiplier * sph_tmp.colour[0];
                                    pixelColour[1] += lightingMultiplier * sph_tmp.colour[1];
                                    pixelColour[2] += lightingMultiplier * sph_tmp.colour[2];
                                }
                                else if (plane_tmp.dist != -1)
                                {
                                    reflectionValue = reflectionValue * plane_tmp.attributes[0];
                                    lightingMultiplier = scene.light.intensity * reflectionValue * 0.65f * intensity;
                                    pixelColour[0] += lightingMultiplier * plane_tmp.colour[0];
                                    pixelColour[1] += lightingMultiplier * plane_tmp.colour[1];
                                    pixelColour[2] += lightingMultiplier * plane_tmp.colour[2];
                                }

                            }
                            temp_pixel_val = true;

                            reflectDepth++;
                        }
                        else
                            reflectDepth = 99;
                    }

                    //inversion from crossing plane
                    //pixelb[x_dir + 1000, -y_dir + 1000, 0] = Convert.ToByte(pixelColour[2]);
                    //pixelb[x_dir + 1000, -y_dir + 1000, 1] = Convert.ToByte(pixelColour[1]);
                    //pixelb[x_dir + 1000, -y_dir + 1000, 2] = Convert.ToByte(pixelColour[0]);
                    pixelb[-y_dir + 1000, x_dir + 1000 , 0] = Convert.ToByte(pixelColour[2]);
                    pixelb[-y_dir + 1000, x_dir + 1000, 1] = Convert.ToByte(pixelColour[1]);
                    pixelb[-y_dir + 1000, x_dir + 1000, 2] = Convert.ToByte(pixelColour[0]);
                }
            }
        }


        /*
        public static void setBmp()
        {


            PixelFormat pxf = PixelFormat.Format32bppRgb;
            // Lock the bitmap's bits.
            Rectangle rect = new Rectangle(0, 0, bmpPic.Width, bmpPic.Height);
            BitmapData bmpData = bmpPic.LockBits(rect, ImageLockMode.WriteOnly, pxf);

            //first address
            IntPtr ptr = bmpData.Scan0;

            // Declare an array to hold the bytes of the bitmap. 
            // int numBytes = bmp.Width * bmp.Height * 3; 
            int numBytes = bmpPic.Width * 4;
            byte[] rgbValues = new byte[numBytes];


            for (int j = 0; j <= 2000; j++)
            {

                for (int i = 0; i <= 2000; i++)
                {


                        rgbValues[ i * 4 + 2] = Convert.ToByte(pixel[i, j, 0]);
                        rgbValues[i * 4 + 1] = Convert.ToByte(pixel[i, j, 1]);
                        rgbValues[i * 4] = Convert.ToByte(pixel[i, j, 2]);

                }

                // Copy the RGB values back to the bitmap
                Marshal.Copy(rgbValues, 0, ptr+8004*j, numBytes);
            }
            
            // Unlock the bits.
            bmpPic.UnlockBits(bmpData);

        }
        */
    }
}

