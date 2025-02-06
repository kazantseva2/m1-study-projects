import org.hw3.PageVisitCounter;
import org.junit.Test;

import java.util.*;

import static org.hamcrest.CoreMatchers.hasItems;
import static org.junit.Assert.*;

public class PageVisitCounterTest {

    @Test
    public void testCorrectness() {
        List<String> pages = Arrays.asList("page1", "page2", "page3", "page4", "page5");
        PageVisitCounter counter = new PageVisitCounter(pages);
        int pageRequestCount = 10;

        counter.run(pageRequestCount);

        // Проверка числа зарегистрированных страниц
        int pagesNumber = counter.pageVisitMap.keySet().size();
        assertEquals(pagesNumber, 5);

        // Проверка имен зарегистрированных страниц
        for (String page : pages) {
            assertThat(counter.pageVisitMap.keySet(), hasItems(page));
        }

        // Проверка, что для каждой страницы количество запросов корректно
        for (String page : pages) {
            assertEquals("Incorrect count for page: " + page, pageRequestCount,
                    counter.pageVisitMap.get(page).get());
        }

        // Проверка, что все потоки завершили выполнение
        for (Thread visitor : counter.visitors) {
            assertTrue("Thread did not complete: " + visitor.getName(),
                    visitor.getState() == Thread.State.TERMINATED);
        }

        // Проверка количества потоков в списке
        assertEquals("Unexpected number of threads",
                pages.size() * pageRequestCount, counter.visitors.size());

    }

    @Test
    public void testPerformance() {
        int pageRequestCount = 1000;
        int numberOfPages = 100;

        // Генерируем список страниц
        List<String> pages = new ArrayList<>();
        for (int i = 1; i <= numberOfPages; i++) {
            pages.add("page" + i);
        }

        PageVisitCounter counter = new PageVisitCounter(pages);


        long startTime = System.nanoTime();

        counter.run(pageRequestCount);

        long endTime = System.nanoTime();

        long duration = endTime - startTime;

        // Время выполнения не должно быть больше 7 секунд
        System.out.println("Test executed in: " + duration / 1_000_000 + " ms");
        long time = 1000000 * 7;
        System.out.println("Max time: " + time / 1000000 + " seconds");
        assertTrue("Test took too long", duration / 1000 < time);
    }
}
