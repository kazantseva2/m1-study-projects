package org.vmk.dep508.threads;

/**
 * Created by Asus on 10/1/2018.
 */
public class YieldExample {
    public static void main(String[] args) {
        Thread t = new MyThread();
        t.start();

        for(int i = 0;i < 5; i++) {
            Thread.yield();
            System.out.println("Control in "+Thread.currentThread().getName());
        }
    }
}

class MyThread extends Thread {

    @Override
    public void run() {
        for(int i = 0;i < 5; i++) {
            System.out.println("Control in "+Thread.currentThread().getName());
        }
    }
}
