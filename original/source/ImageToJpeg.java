import com.sun.image.codec.jpeg.*;
import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.awt.image.*;
import java.io.*;

/**
 * A class that is used to convert an image into jpeg format.  This code was derived from a developer.com article by Benoit Marchal "JDK 1.2 does JPEG" 11/25/1998
 * @author Benoit Marchal
 * @author Sriram Kosuri
 * @version 1.0
 */
public class ImageToJpeg{
	
	/*
	 * A static method to encode an image into high quality jpeg format (no compression) into an output stream
	 * @param img The BufferedImage that is to be converted to jpeg
	 * @param out The output stream which the jpeg will be output to.  Usually a FileOutputStream.
	 */
    public static void EncodeIt(BufferedImage img, OutputStream out) throws IOException {
        JPEGImageEncoder encoder = JPEGCodec.createJPEGEncoder(out);
        JPEGEncodeParam param = encoder.getDefaultJPEGEncodeParam(img);
        param.setQuality(1.0f,true);
        encoder.encode(img,param);
    }       
}