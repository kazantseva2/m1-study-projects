package org.vmk.dep508.threads;

import java.util.concurrent.*;

/**
 * Created by Asus on 10/1/2018.
 */
public class ThreadPoolExecutor {
    public static void main(String[] args) throws ExecutionException, InterruptedException {
        ExecutorService service = Executors.newFixedThreadPool(10);

        Future<Integer>[] results = new Future[30];
        for(int i = 0; i < 30; i++) {
            IntegerCallable c = new IntegerCallable(i);
            results[i] = service.submit(c);
        }

        for(int i = 0; i < results.length; i++) {
            Future<Integer> f = (Future<Integer>)results[i];
            System.out.println("Future result in "+i+" is "+f.get());
        }
    }
}

class IntegerCallable implements Callable<Integer> {

    Integer i;

    public IntegerCallable(Integer i) {
        this.i = i;
    }

    public Integer call() throws Exception {
        System.out.println("Within callable - "+i);
        return ++i;
    }
}