package org.hw3;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;


public class PageVisitCounter {
    public final List<VisitorThread> visitors = new ArrayList<>();
    public final Map<String, AtomicInteger> pageVisitMap = new ConcurrentHashMap<>();

    public PageVisitCounter(List<String> pages) {
        for (String page : pages) {
            this.pageVisitMap.put(page, new AtomicInteger(0));
        }
    }

    public void run(int pageRequestCount) {
        for (String page : pageVisitMap.keySet()) {
            for (int i = 0; i < pageRequestCount; i++) {
                visitors.add(new VisitorThread(page, pageVisitMap));
            }
        }

        for (VisitorThread visitor: visitors){
            visitor.start();
        }

        for (Thread visitor : visitors) {
            try {
                visitor.join();
            } catch(InterruptedException e){}
        }
    }
}