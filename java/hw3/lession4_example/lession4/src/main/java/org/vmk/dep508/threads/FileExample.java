package org.vmk.dep508.threads;

import java.io.File;

/**
 * Created by Asus on 10/1/2018.
 */
public class FileExample {
    public static void main(String[] args) {
        File f = new File("D:\\tmp\\");
        System.out.println(f.exists());
        System.out.println(f.isDirectory());
        for(File file : f.listFiles()){
            System.out.println(file.getName());
        }

    }
}
