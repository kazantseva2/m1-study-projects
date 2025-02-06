package org.hw3;

import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

class VisitorThread extends Thread {
    private final String pageName;
    private final Map<String, AtomicInteger> pageVisitMap;

    VisitorThread(String pageName, Map<String, AtomicInteger> pageVisitMap) {
        this.pageName = pageName;
        this.pageVisitMap = pageVisitMap;
    }

    @Override
    public void run() {
        pageVisitMap.get(pageName).incrementAndGet();
    }
}