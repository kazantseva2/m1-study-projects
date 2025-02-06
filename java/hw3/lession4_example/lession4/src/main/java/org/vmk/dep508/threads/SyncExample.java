package org.vmk.dep508.threads;

/**
 * Created by Asus on 10/1/2018.
 */
public class SyncExample {
    public static void main(String[] args) {
        Incrementor inc = new Incrementor();

        IncrementThread t1 = new IncrementThread();
        IncrementThread t2 = new IncrementThread();
        IncrementThread t3 = new IncrementThread();
        t1.setInc(inc);
        t2.setInc(inc);
        t3.setInc(inc);


        t1.start();
        t2.start();
        t3.start();

        synchronized (inc) {
            inc.notifyAll();
        }

        try {
            t1.join();
            t2.join();
            t3.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }


        System.out.println(inc.getI());
    }
}

class IncrementThread extends Thread {

    Incrementor inc;

    public void setInc(Incrementor inc) {
        this.inc = inc;
    }

    @Override
    public void run() {
        for(int i = 0; i < 1000; i++) {
            //inc.increment();
            //inc.syncIncrement();
            inc.incrementWithWait();
        }
    }
}

class Incrementor {
    private int i = 0;
    private Object o = new Object();

    public int getI() {
        return i;
    }

    public void increment(){
        i++;

    }

    public synchronized void syncIncrement(){
        i++;
    }

    public synchronized void incrementWithWait () {
        try {
            wait(1000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        i++;
        notifyAll();
    }
}
