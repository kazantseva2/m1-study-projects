package org.vmk.dep508.collections.auction;

import org.junit.*;

import java.math.BigDecimal;
import java.util.List;

import static org.junit.Assert.*;
import static org.hamcrest.CoreMatchers.*;

public class AuctionImplTest {

    Auction auction;
    String product4thEdition = "Thinking Java 4th edition";
    String product3rdEdition = "Thinking Java 3rd edition";
    BigDecimal defaultPrice = BigDecimal.TEN;

    /*Вызыватется при инициализации класса AuctionImplTest*/
    @BeforeClass
    public static void setupClass(){

    }

    /*Вызыватется перед вызовом каждого метода помеченного аннотацией @Test*/
    @Before
    public void setup(){
        auction = new AuctionImpl();
        auction.placeProduct(product3rdEdition, defaultPrice);
        auction.placeProduct(product4thEdition, defaultPrice);
    }

    /*Вызыватется после вызова каждого метода помеченного аннотацией @Test*/
    @After
    public void clear() {
        auction = null;
    }

    /*Вызывается после вызова всех тестовых методов*/
    @AfterClass
    public static void releaseRecources() {

    }

    @Test(expected = ProductNotFoundException.class)
    public void placeProduct() throws Exception {
        List<String> products = auction.getProducts();

        assertThat(products, hasItems(product3rdEdition, product4thEdition));
        assertEquals(products.size(), 2);

        String product5thEdition = "Thinking Java 5th edition";
        assertThat(products, not(hasItem(product5thEdition)));

        assertEquals(auction.getProductPrice(product3rdEdition), defaultPrice);

        BigDecimal notExistingproductPrice = auction.getProductPrice(product5thEdition);
    }

    @Test(expected = ProductNotFoundException.class)
    public void addBid() throws Exception {
        String user1 = "user1";
        String user2 = "user2";
        BigDecimal price1 =  new BigDecimal(100);
        BigDecimal price2 =  new BigDecimal(1);

        auction.addBid(user1, product3rdEdition, price1);
        auction.addBid(user2, product3rdEdition, price2);
        assertEquals(auction.getProductPrice(product3rdEdition), price1);

        auction.addBid(user2, product4thEdition, price2);
        assertEquals(auction.getProductPrice(product4thEdition), defaultPrice);

        String product5thEdition = "Thinking Java 5th edition";
        auction.addBid(user1, product5thEdition, price1);
    }

    @Test
    public void removeBid() throws Exception {
        String user1 = "user1";
        String user2 = "user2";
        BigDecimal price1 =  new BigDecimal(100);
        BigDecimal price2 =  new BigDecimal(50);


        auction.addBid(user1, product3rdEdition, price1);
        assertEquals(auction.getProductPrice(product3rdEdition), price1);

        auction.removeBid(user2, product3rdEdition);
        assertEquals(auction.getProductPrice(product3rdEdition), price1);

        auction.removeBid(user1, product4thEdition);
        assertEquals(auction.getProductPrice(product3rdEdition), price1);

        auction.addBid(user2, product3rdEdition, price2);

        auction.removeBid(user1, product3rdEdition);
        assertEquals(auction.getProductPrice(product3rdEdition), price2);

        auction.removeBid(user2, product3rdEdition);
        assertEquals(auction.getProductPrice(product3rdEdition), defaultPrice);
    }

    @Test(expected = ProductNotFoundException.class)
    public void sellProduct() throws Exception {
        String user1 = "user1";
        BigDecimal price1 =  new BigDecimal(1);
        BigDecimal price2 =  new BigDecimal(100);

        auction.addBid(user1, product3rdEdition, price1);
        assertFalse(auction.sellProduct(product3rdEdition));

        auction.addBid(user1, product3rdEdition, price2);
        assertTrue(auction.sellProduct(product3rdEdition));

        String product5thEdition = "Thinking Java 5th edition";
        boolean resSell = auction.sellProduct(product5thEdition);
    }

}
