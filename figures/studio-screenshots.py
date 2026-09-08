import asyncio, time
from playwright.async_api import async_playwright
URL="https://users.ox.ac.uk/~coml1155/mosaics-studio/"
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch()
        ctx=await b.new_context(viewport={"width":1600,"height":1000}, device_scale_factor=2)
        page=await ctx.new_page()
        page.on("console", lambda m: None)
        await page.goto(URL, wait_until="load")
        await page.wait_for_timeout(15000)
        await page.screenshot(path="shots/step1-start.png")
        await page.screenshot(path="shots/step1-start-full.png", full_page=True)
        await page.get_by_role("button", name="Build the structure").click()
        await page.wait_for_timeout(10000)
        await page.screenshot(path="shots/step2-prepare.png")
        await page.screenshot(path="shots/step2-prepare-full.png", full_page=True)
        # protocol choices on the run page
        await page.get_by_role("button", name="Set up the run", exact=True).click()
        await page.wait_for_timeout(3000)
        for s in await page.locator("select").all():
            try:
                if await s.is_visible():
                    print("SELECT visible:", [o[:70] for o in await s.locator('option').all_inner_texts()][:20])
            except Exception: pass
        await page.screenshot(path="shots/step3-run.png")
        await page.screenshot(path="shots/step3-run-full.png", full_page=True)
        await page.get_by_role("button", name="Run the engine").click()
        t0=time.time()
        # wait for completion: the 'Run again' button becomes visible
        for i in range(120):
            await page.wait_for_timeout(2000)
            try:
                if await page.get_by_role("button", name="Run again").is_visible():
                    break
            except Exception: pass
        print("run finished after", round(time.time()-t0,1), "s")
        await page.wait_for_timeout(3000)
        await page.screenshot(path="shots/step4-understand.png")
        await page.screenshot(path="shots/step4-understand-full.png", full_page=True)
        txt=await page.locator("body").inner_text()
        open("shots/step4-text.txt","w").write(txt)
        await b.close()
asyncio.run(main())
