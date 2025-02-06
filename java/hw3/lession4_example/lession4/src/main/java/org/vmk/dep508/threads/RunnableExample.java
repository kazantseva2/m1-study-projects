package org.vmk.dep508.threads;

/**
 * Hello world!
 *
 */
public class RunnableExample
{

    public static void main( String[] args )
    {
        final Integer i = new Integer(1);
        Thread th = new Thread(new Runnable() {
            public void run() {
                System.out.println("Hello from thread "+Thread.currentThread().getName());

                //System.out.println(i);

                //thread is sleeping for 1 sec
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }

                System.out.println("Thread "+Thread.currentThread().getName()+" is completing it's work");
            }
        });

        th.setName("Our test thread");
        th.start();

        //witing while thread is working
        try {
            th.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        System.out.println("Thread "+th.getName()+" completed it's work");
    }

}
